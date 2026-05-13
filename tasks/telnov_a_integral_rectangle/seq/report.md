# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) - SEQ

- Student: Тельнов А., group 3823Б1ФИ1
- Technology: SEQ
- Variant: 9

## 1. Introduction

Последовательная версия используется как базовая линия производительности. Она нужна для проверки математической корректности метода прямоугольников и для дальнейшего расчёта ускорения параллельных реализаций.

## 2. Problem Statement

Вход: `InType = std::pair<int, int>`, где `first` — число разбиений `n`, `second` — размерность `dimensions`.

Выход: `OutType = double`, приближённое значение интеграла функции `f(x1, ..., xd) = x1 + ... + xd` на области `[0, 1]^d`.

Ограничения: `n > 0`, `dimensions > 0`.

## 3. Baseline Algorithm (Sequential)

SEQ-версия выполняет полный перебор всех `n^dimensions` точек в одном потоке. Линейный индекс `idx` преобразуется в координаты многомерной сетки через остаток и деление на `n`. Для каждой координаты используется центр ячейки: `x = (coord_index + 0.5) * h`.

Значение функции в точке равно сумме координат. После обхода всех точек накопленная сумма умножается на объём ячейки `h^dimensions`.

Эта версия не использует потоки, MPI, TBB или OpenMP. Поэтому время SEQ считается baseline для расчёта `speedup` во всех остальных отчётах.

## 4. Parallelization Scheme

Параллелизм отсутствует. Все операции выполняются последовательно в одном потоке. `workers = 1`.

## 5. Implementation Details

- Файлы: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`.
- Класс: `TelnovAIntegralRectangleSEQ`.
- `ValidationImpl()` проверяет положительность `n` и `dimensions`.
- `PreProcessingImpl()` обнуляет выходное значение.
- `RunImpl()` выполняет основной расчёт.
- `PostProcessingImpl()` завершает задачу без дополнительных преобразований.

Память под многомерную сетку не выделяется. Координаты вычисляются на лету из линейного индекса, поэтому дополнительная память имеет порядок `O(1)`.

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
bool TelnovAIntegralRectangleSEQ::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / n;

  auto total_points = static_cast<int64_t>(std::pow(n, d));

  double result = 0.0;

  for (int64_t idx = 0; idx < total_points; idx++) {
    int64_t tmp = idx;
    double f_value = 0.0;

    for (int dim = 0; dim < d; dim++) {
      int coord_index = static_cast<int>(tmp % n);
      tmp /= n;

      double x = a + ((coord_index + 0.5) * h);
      f_value += x;
    }

    result += f_value;
  }

  result *= std::pow(h, d);

  GetOutput() = result;
  return true;
}
```
