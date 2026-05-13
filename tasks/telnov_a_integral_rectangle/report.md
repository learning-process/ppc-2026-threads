# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников)

- Student: Тельнов А., group 3823Б1ФИ1
- Technology: SEQ, OMP, TBB, STL, ALL
- Variant: 9

## 1. Introduction

Работа посвящена численному вычислению многомерного интеграла методом прямоугольников. Такой подход применяется, когда аналитическое интегрирование неудобно или невозможно, а область интегрирования можно представить как сетку. В работе реализованы пять вариантов одной задачи: последовательный baseline и параллельные версии на OMP, TBB, STL и ALL.

Ожидаемый результат работы — корректное приближённое значение интеграла и сравнение времени выполнения разных технологий параллелизма на одинаковой постановке задачи.

## 2. Problem Statement

Входные данные задаются типом `InType = std::pair<int, int>`:

- `n` — количество прямоугольников по каждой координате;
- `dimensions` — размерность интеграла.

Выходные данные: `OutType = double`, приближённое значение интеграла.

Область интегрирования — единичный гиперкуб:

```txt
[0, 1]^dimensions
```

Интегрируемая функция:

```txt
f(x1, x2, ..., xd) = x1 + x2 + ... + xd
```

Метод использует центры элементарных ячеек. Шаг сетки `h = 1 / n`, общее число точек `N = n^dimensions`. Ограничения: `n > 0`, `dimensions > 0`.

## 3. Baseline Algorithm (Sequential)

Последовательная версия служит baseline для всех остальных реализаций. Алгоритм перебирает линейные индексы от `0` до `n^dimensions - 1`. Каждый линейный индекс переводится в набор координат с помощью деления и остатка по `n`. Для каждой координаты берётся центр соответствующей ячейки: `x = (coord_index + 0.5) * h`.

Затем вычисляется значение функции как сумма координат. Все значения суммируются, после чего итоговая сумма умножается на объём элементарной ячейки: `result = sum * h^dimensions`.

## 4. Parallelization Scheme

- SEQ: один поток, параллелизм не используется.
- OMP: параллельный цикл по линейному диапазону индексов с редукцией суммы.
- TBB: `oneapi::tbb::parallel_reduce` по `blocked_range<int64_t>`.
- STL: ручное разбиение диапазона между `std::thread`, частичные суммы и последующее объединение после `join()`.
- ALL: MPI делит общий диапазон между rank-ами, внутри каждого rank-а используются STL-потоки, TBB для подсчёта поддиапазонов и OpenMP для редукции частичных сумм; глобальное объединение выполняется через `MPI_Allreduce`.

## 5. Implementation Details

Код задачи расположен в папке `tasks/telnov_a_integral_rectangle`:

- `common/include/common.hpp` — типы `InType`, `OutType`, `TestType`, `BaseTask`;
- `seq/src/ops_seq.cpp` — последовательный baseline;
- `omp/src/ops_omp.cpp` — OpenMP-версия;
- `tbb/src/ops_tbb.cpp` — TBB-версия;
- `stl/src/ops_stl.cpp` — STL-версия;
- `all/src/ops_all.cpp` — гибридная версия ALL;
- `tests/functional/main.cpp` — функциональные тесты;
- `tests/performance/main.cpp` — performance-тесты.

Память используется экономно: основная реализация не хранит всю многомерную сетку, а вычисляет координаты точки из линейного индекса на лету.

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
| tbb | 2 | 0.122722 | 1.53 | 76.37% |
| tbb | 4 | 0.062125 | 3.02 | 75.43% |
| tbb | 8 | 0.044217 | 4.24 | 52.99% |
| stl | 2 | 0.098049 | 1.91 | 95.58% |
| stl | 4 | 0.070281 | 2.67 | 66.67% |
| stl | 8 | 0.054564 | 3.44 | 42.94% |
| all | 2 × 1 | 0.114299 | 1.64 | 82.00% |
| all | 2 × 2 | 0.068596 | 2.73 | 68.31% |
| all | 4 × 2 | 0.067989 | 2.76 | 34.46% |

## 8. Conclusions

Реализация вычисляет многомерный интеграл методом прямоугольников без хранения полной сетки в памяти. Последовательная версия используется как baseline, а параллельные версии сравниваются с ней по времени выполнения, ускорению и эффективности. Основное ограничение масштабируемости связано с ростом числа точек `n^dimensions`, накладными расходами на управление параллелизмом и редукцией частичных сумм.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

Ниже приведены основные фрагменты `RunImpl()` для всех реализаций.

### SEQ RunImpl

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

### OMP RunImpl

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

### TBB RunImpl

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

### STL RunImpl

```cpp
bool TelnovAIntegralRectangleSTL::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  int thread_count = ppc::util::GetNumThreads();
  thread_count = std::max(1, std::min(thread_count, static_cast<int>(total_points)));

  std::vector<std::thread> threads(thread_count);
  std::vector<double> partial_sums(thread_count, 0.0);

  const int64_t block = total_points / thread_count;
  const int64_t remainder = total_points % thread_count;

  auto calculate_part = [n, d, a, h](int64_t begin, int64_t end) {
    double local_sum = 0.0;

    for (int64_t idx = begin; idx < end; ++idx) {
      int64_t current = idx;
      double f_value = 0.0;

      for (int dim = 0; dim < d; ++dim) {
        const int coord_index = static_cast<int>(current % n);
        current /= n;

        const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
        f_value += x;
      }

      local_sum += f_value;
    }

    return local_sum;
  };

  int64_t begin = 0;
  for (int i = 0; i < thread_count; ++i) {
    const int64_t current_block = block + (i < remainder ? 1 : 0);
    const int64_t end = begin + current_block;

    threads[i] = std::thread([&, i, begin, end]() { partial_sums[i] = calculate_part(begin, end); });

    begin = end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  double result = 0.0;
  for (const auto &value : partial_sums) {
    result += value;
  }

  GetOutput() = result * std::pow(h, d);
  return true;
}
```

### ALL RunImpl

```cpp
bool TelnovAIntegralRectangleALL::RunImpl() {
  const int n = GetInput().first;
  const int dimensions = GetInput().second;

  const double h = 1.0 / static_cast<double>(n);
  const int64_t total_points = PowInt(n, dimensions);

  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int64_t base_block = total_points / size;
  const int64_t remainder = total_points % size;

  const int64_t rank_begin = (static_cast<int64_t>(rank) * base_block) + std::min<int64_t>(rank, remainder);
  const int64_t rank_size = base_block + (rank < remainder ? 1 : 0);

  int thread_count = ppc::util::GetNumThreads();
  thread_count = std::max(1, std::min(thread_count, static_cast<int>(std::max<int64_t>(1, rank_size))));

  std::vector<std::thread> threads(thread_count);
  std::vector<double> thread_sums(thread_count, 0.0);

  const int64_t thread_block = rank_size / thread_count;
  const int64_t thread_remainder = rank_size % thread_count;

  int64_t current_begin = rank_begin;

  for (int thread_id = 0; thread_id < thread_count; ++thread_id) {
    const int64_t current_block = thread_block + (thread_id < thread_remainder ? 1 : 0);
    const int64_t current_end = current_begin + current_block;

    threads[thread_id] = std::thread([&, thread_id, current_begin, current_end]() {
      thread_sums[thread_id] = CalculateRange(current_begin, current_end, n, dimensions, h);
    });

    current_begin = current_end;
  }

  for (auto &thread : threads) {
    thread.join();
  }

  double local_sum = 0.0;

#pragma omp parallel for default(none) shared(thread_sums, thread_count) reduction(+ : local_sum) \
    num_threads(thread_count)
  for (int i = 0; i < thread_count; ++i) {
    local_sum += thread_sums[i];
  }

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  GetOutput() = global_sum;

  for (int i = 0; i < dimensions; ++i) {
    GetOutput() *= h;
  }

  return true;
}
```
