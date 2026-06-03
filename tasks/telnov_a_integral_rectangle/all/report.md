# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) - ALL

- Student: Тельнов А., group 3823Б1ФИ1
- Technology: ALL
- Variant: 9

## 1. Introduction

ALL-версия объединяет несколько технологий параллелизма. Цель этой реализации — распределить общий диапазон
точек между MPI-процессами, а внутри каждого процесса дополнительно использовать потоки и локальные редукции.

## 2. Problem Statement

Задача совпадает с остальными версиями: вычислить интеграл функции `f(x1, ..., xd) = x1 + ... + xd` на `[0,
1]^d` методом прямоугольников.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Он выполняет один последовательный обход всех `n^dimensions` точек. Для
ALL-версии этот результат используется как эталон для проверки корректности и расчёта ускорения.

## 4. Parallelization Scheme

ALL-реализация использует гибридную схему:

1. **MPI:** общий диапазон `[0, total_points)` делится между rank-ами. Каждый rank получает непрерывный блок
индексов.
2. **STL:** внутри rank-а локальный блок дополнительно делится между `std::thread`.
3. **TBB:** каждый STL-поток вызывает `CalculateRange`, где поддиапазон обрабатывается через
`oneapi::tbb::parallel_reduce` и `blocked_range<int64_t>`.
4. **OpenMP:** частичные суммы STL-потоков объединяются через `#pragma omp parallel for ... reduction(+ :
local_sum)`.
5. **MPI-синхронизация:** `MPI_Allreduce` складывает `local_sum` всех rank-ов и возвращает одинаковую
`global_sum` всем процессам.

Конфигурация задаётся как `workers = ranks × threads`. Например, при запуске `mpiexec -n 4` и
`PPC_NUM_THREADS=2` конфигурация равна `4 × 2`, а `workers = 8`.

Смысл `MPI_Allreduce` — синхронизировать итоговую сумму между процессами. После этой операции каждый rank
имеет одно и то же глобальное значение интегральной суммы.

## 5. Implementation Details

- Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.
- Класс: `TelnovAIntegralRectangleALL`.
- `PowInt` вычисляет `n^dimensions` без хранения сетки.
- `CalculatePointValue` вычисляет значение функции в одной точке.
- `CalculateRange` выполняет TBB-редукцию на поддиапазоне.
- `RunImpl` отвечает за MPI-разбиение, запуск STL-потоков, OpenMP-редукцию и финальный `MPI_Allreduce`.

Память: на каждом rank-е хранится массив STL-потоков и массив `thread_sums`. Полная сетка не создаётся.

Corner cases:

- Если `total_points` не делится на число rank-ов, остаток распределяется между первыми rank-ами.
- Если точек меньше, чем потоков, число потоков ограничивается размером локального блока.

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
| all | 2 × 1 | 0.114299 | 1.64 | 82.00% |
| all | 2 × 2 | 0.068596 | 2.73 | 68.31% |
| all | 4 × 2 | 0.067989 | 2.76 | 34.46% |

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
