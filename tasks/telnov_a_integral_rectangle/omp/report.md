# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) - OMP

- Student: Тельнов А., group 3823Б1ФИ1
- Technology: OMP
- Variant: 9

## 1. Introduction

OpenMP-версия реализует параллельный перебор независимых точек сетки. Метод прямоугольников хорошо подходит для такой схемы, потому что значение функции в каждой точке не зависит от других точек.

## 2. Problem Statement

Задача совпадает с SEQ-версией: вычислить интеграл функции `f(x1, ..., xd) = x1 + ... + xd` на `[0, 1]^d` при заданных `n` и `dimensions`.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Он выполняет последовательный обход `n^dimensions` точек, вычисляет сумму координат в центре каждой ячейки и умножает накопленную сумму на `h^dimensions`.

## 4. Parallelization Scheme

В OMP-версии параллелится внешний цикл по линейному диапазону индексов:

```cpp
#pragma omp parallel for default(none) reduction(+ : result) shared(total_points, n, d, a, h)
```

- `shared`: `total_points`, `n`, `d`, `a`, `h` — эти значения только читаются потоками.
- `private`: `idx`, `tmp`, `f_value`, `dim`, `coord_index`, `x` — локальные переменные итераций цикла.
- `reduction(+ : result)`: каждый поток накапливает свою локальную сумму, затем OpenMP объединяет частичные суммы.
- `schedule`: явно не указан, поэтому используется расписание по умолчанию компилятора OpenMP. Для данной задачи это допустимо, потому что стоимость обработки каждой точки примерно одинакова.

Синхронизация выполняется OpenMP в конце параллельного цикла при объединении редукции.

## 5. Implementation Details

- Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.
- Класс: `TelnovAIntegralRectangleOMP`.
- Основное отличие от SEQ — директива `parallel for` с редукцией суммы.
- Гонок данных нет, так как общий результат обновляется через `reduction`, а остальные переменные внутри цикла являются локальными.

Память: дополнительно используется только внутреннее хранилище OpenMP для частичных сумм редукции. Полная сетка в памяти не создаётся.

## 6. Experimental Setup

**Аппаратное обеспечение:**

- **CPU:** 12th Gen Intel(R) Core(TM) i5-12450H (2.00 GHz, 8 ядер / 12 потоков)
- **RAM:** 16 ГБ
- **OS:** Windows 11 Pro x64
- **MPI:** Microsoft MPI (MS-MPI) 10.1

**Инструменты:**

- **Сборка:** CMake
- **Компилятор:** MSVC 19.x
- **Конфигурация:** Release

**Окружение:**

- **PPC_NUM_THREADS:** задаёт число потоков для OMP, TBB, STL и потоковой части ALL.
- **PPC_NUM_PROC / mpiexec -n:** задаёт число MPI-процессов для ALL.
- Для ALL конфигурация записывается в формате `ranks × threads`.

**Генерация данных:**

- Тесты генерируют входные данные автоматически.
- Для performance-теста используется вход `InType{50, 4}`: 50 разбиений по каждой координате и размерность 4.
- Внешние файлы с данными не используются.

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами из `tests/functional/main.cpp`. Для набора параметров `(n, dimensions)` вычисленное значение сравнивается с аналитическим результатом:

```txt
I = dimensions / 2
```

Это значение получается для функции `f(x1, ..., xd) = x1 + ... + xd` на единичном гиперкубе `[0, 1]^d`: интеграл каждой координаты равен `1/2`, поэтому сумма по `d` координатам равна `d / 2`. В тестах используется допустимое отклонение, зависящее от числа разбиений `n`.

### 7.2 Performance

Используемые обозначения:

```txt
time — время выполнения performance-теста;
speedup = time_seq / time_mode;
efficiency = speedup / workers;
workers — количество исполнителей: потоков для OMP/TBB/STL, ranks × threads для ALL.
```
| Mode | Count | Time, s | Speedup | Efficiency |
|-------------|-------|---------|---------|------------|
| seq | 1 | 0.187439 | 1.00 | N/A |
| omp | 2 | 0.047289 | 3.96 | 198.18% |
| omp | 4 | 0.051861 | 3.61 | 90.36% |
| omp | 8 | 0.044240 | 4.24 | 52.96% |

В таблицу необходимо подставить реальные значения времени из локального запуска `ppc_perf_tests`. После этого `speedup` и `efficiency` рассчитываются по формулам выше. Без фактических замеров не делается утверждение об оптимальности реализации.

## 8. Conclusions

Реализация вычисляет многомерный интеграл методом прямоугольников без хранения полной сетки в памяти. Последовательная версия используется как baseline, а параллельные версии сравниваются с ней по времени выполнения, ускорению и эффективности. Основное ограничение масштабируемости связано с ростом числа точек `n^dimensions`, накладными расходами на управление параллелизмом и редукцией частичных сумм.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

Основной фрагмент `RunImpl()`:

```cpp
bool TelnovAIntegralRectangleOMP::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  double result = 0.0;

#pragma omp parallel for default(none) reduction(+ : result) shared(total_points, n, d, a, h)
  for (int64_t idx = 0; idx < total_points; idx++) {
    int64_t tmp = idx;
    double f_value = 0.0;

    for (int dim = 0; dim < d; dim++) {
      const int coord_index = static_cast<int>(tmp % n);
      tmp /= n;

      const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
      f_value += x;
    }

    result += f_value;
  }

  result *= std::pow(h, d);

  GetOutput() = result;
  return true;
}
```
