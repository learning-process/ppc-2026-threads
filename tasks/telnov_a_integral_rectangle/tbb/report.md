# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) - TBB

- Student: Тельнов А., group 3823Б1ФИ1
- Technology: TBB
- Variant: 9

## 1. Introduction

TBB-версия использует задачно-ориентированный подход oneAPI Threading Building Blocks. Вместо явного создания
потоков общий диапазон индексов передаётся планировщику TBB, который разбивает его на блоки и распределяет
работу между потоками.

## 2. Problem Statement

Требуется вычислить интеграл функции `f(x1, ..., xd) = x1 + ... + xd` на единичном гиперкубе `[0, 1]^d`
методом прямоугольников.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`: один поток обходит `n^dimensions` точек, переводит линейный индекс в
координаты, суммирует значения функции и умножает результат на объём ячейки.

## 4. Parallelization Scheme

Используется `oneapi::tbb::parallel_reduce`:

```cpp
oneapi::tbb::parallel_reduce(
    oneapi::tbb::blocked_range<int64_t>(0, total_points),
    0.0,
    body,
    join);
```

- `blocked_range<int64_t>` задаёт полуинтервал линейных индексов `[begin, end)`.
- `grainsize` явно не задан в коде, поэтому используется значение по умолчанию для `blocked_range`.
- `partitioner` явно не передаётся; применяется стандартное поведение TBB для выбранной перегрузки
  `parallel_reduce`.
- Каждый поддиапазон возвращает локальную сумму.
- Объединение локальных сумм выполняется функцией `join`, которая складывает два значения `double`.

Контроль конкуренции выполняется через окружение и настройки тестового фреймворка PPC (`PPC_NUM_THREADS`). В
самой задаче `global_control` не создаётся, поэтому код не меняет глобальные настройки TBB для других тестов.

## 5. Implementation Details

- Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.
- Класс: `TelnovAIntegralRectangleTBB`.
- В теле `parallel_reduce` каждая задача обрабатывает свой `blocked_range` и вычисляет частичную сумму.
- Гонок данных нет, так как `local_sum` принадлежит конкретной задаче, а общий результат формируется только
  через редукцию TBB.

Память: не создаётся массив точек сетки; TBB хранит только служебные задачи и локальные суммы.

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

Корректность проверялась функциональными тестами из `tests/functional/main.cpp`. Для набора параметров `(n,
dimensions)` вычисленное значение сравнивается с аналитическим результатом:

```txt
I = dimensions / 2
```

Это значение получается для функции `f(x1, ..., xd) = x1 + ... + xd` на единичном гиперкубе `[0, 1]^d`:
интеграл каждой координаты равен `1/2`, поэтому сумма по `d` координатам равна `d / 2`. В тестах используется
допустимое отклонение, зависящее от числа разбиений `n`.

### 7.2 Performance

Используемые обозначения:

```txt
time — время выполнения performance-теста;
speedup = time_seq / time_mode;
efficiency = speedup / workers;
workers — количество исполнителей: потоков для OMP/TBB/STL, ranks × threads для ALL.
```

| Mode | Count | Time, s | Speedup | Efficiency |
| --- | --- | --- | --- | --- |
| seq | 1 | 0.187439 | 1.00 | N/A |
| tbb | 2 | 0.122722 | 1.53 | 76.37% |
| tbb | 4 | 0.062125 | 3.02 | 75.43% |
| tbb | 8 | 0.044217 | 4.24 | 52.99% |

## 8. Conclusions

Реализация вычисляет многомерный интеграл методом прямоугольников без хранения полной сетки в памяти.
Последовательная версия используется как baseline, а параллельные версии сравниваются с ней по времени
выполнения, ускорению и эффективности. Основное ограничение масштабируемости связано с ростом числа точек
`n^dimensions`, накладными расходами на управление параллелизмом и редукцией частичных сумм.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

Основной фрагмент `RunImpl()`:

```cpp
bool TelnovAIntegralRectangleTBB::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  const double result = oneapi::tbb::parallel_reduce(
      oneapi::tbb::blocked_range<int64_t>(0, total_points), 0.0,
      [n, d, a, h](const oneapi::tbb::blocked_range<int64_t> &range, double local_sum) {
        for (int64_t idx = range.begin(); idx != range.end(); ++idx) {
          int64_t tmp = idx;
          double f_value = 0.0;

          for (int dim = 0; dim < d; ++dim) {
            const int coord_index = static_cast<int>(tmp % n);
            tmp /= n;

            const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
            f_value += x;
          }

          local_sum += f_value;
        }
        return local_sum;
      },
      [](double lhs, double rhs) { return lhs + rhs; });

  GetOutput() = result * std::pow(h, d);
  return true;
}
```
