# Поразрядная сортировка для вещественных чисел (тип double) с чётно-нечётным слиянием Бэтчера

- Student: Титаев М., group 3823Б1ФИ1
- Technology: SEQ, OMP, TBB, STL, ALL
- Variant: 20

## 1. Introduction

Работа посвящена сортировке массива вещественных чисел типа `double` с помощью поразрядной сортировки
(radix sort) и последующего применения сети чётно-нечётного слияния Бэтчера. Поразрядная сортировка
обеспечивает линейную по числу элементов сложность за счёт обработки чисел по разрядам, а сеть Бэтчера
служит детерминированной сортирующей сетью, удобной для параллельной реализации. В работе реализованы
пять вариантов одной задачи: последовательный baseline и параллельные версии на OMP, TBB, STL и ALL.

Ожидаемый результат работы — корректно отсортированный массив и сравнение времени выполнения разных
технологий параллелизма на одинаковой постановке задачи.

## 2. Problem Statement

Входные данные задаются типом `InType = std::vector<double>` — массив вещественных чисел.

Выходные данные: `OutType = std::vector<double>` — тот же массив, отсортированный по неубыванию.

Ключевая идея — преобразовать каждое число `double` в беззнаковый 64-битный ключ с сохранением порядка,
отсортировать ключи поразрядной сортировкой и при размере, равном степени двойки, применить сеть Бэтчера.
Ограничение: входной массив не должен быть пустым.

## 3. Baseline Algorithm (Sequential)

Последовательная версия служит baseline для всех остальных реализаций. Алгоритм состоит из трёх этапов.

На первом этапе каждое число `double` преобразуется в 64-битный беззнаковый ключ `uint64_t` так, чтобы
лексикографический порядок ключей совпадал с числовым порядком исходных чисел. Для этого у положительных
чисел инвертируется знаковый бит, а у отрицательных инвертируются все биты.

На втором этапе выполняется поразрядная сортировка ключей по основанию 256 (по 8 бит за проход, всего
8 проходов на 64 бита). Каждый проход использует счётную сортировку: подсчёт количества элементов в
корзинах, префиксные суммы и стабильное распределение по временному массиву.

На третьем этапе ключи преобразуются обратно в `double`. Если размер массива является степенью двойки,
дополнительно применяется сеть чётно-нечётного слияния Бэтчера, реализованная как битоническая
сортировка с обменами по индексам, отличающимся на один бит.

## 4. Parallelization Scheme

- SEQ: один поток, параллелизм не используется.
- OMP: параллельное преобразование ключей, подсчёт корзин с локальными гистограммами по потокам и
  параллельные проходы сети Бэтчера через `#pragma omp parallel for`.
- TBB: `oneapi::tbb::parallel_for` по `blocked_range<size_t>` на этапах преобразования и в сети Бэтчера.
- STL: ручное разбиение диапазона между `std::thread`, локальные гистограммы и объединение после
  `join()`.
- ALL: MPI-окружение для запуска под несколькими процессами, внутри каждого процесса OpenMP-параллелизм
  на этапах поразрядной сортировки и сети Бэтчера.

## 5. Implementation Details

Код задачи расположен в папке `tasks/titaev_m_sortirovka_betchera`:

- `common/include/common.hpp` — типы `InType`, `OutType`, `TestType`, `BaseTask`;
- `seq/src/ops_seq.cpp` — последовательный baseline;
- `omp/src/ops_omp.cpp` — OpenMP-версия;
- `tbb/src/ops_tbb.cpp` — TBB-версия;
- `stl/src/ops_stl.cpp` — STL-версия;
- `all/src/ops_all.cpp` — гибридная версия ALL;
- `tests/functional/main.cpp` — функциональные тесты;
- `tests/performance/main.cpp` — performance-тесты.

Память используется экономно: для поразрядной сортировки выделяется один временный массив ключей и массив
счётчиков фиксированного размера на 256 корзин.

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
- Для performance-теста используется массив из `2^20` элементов с убывающими значениями.
- Внешние файлы с данными не используются.

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами из `tests/functional/main.cpp`. Выходной массив
проверяется на полную упорядоченность по неубыванию, а его размер сравнивается с размером входа.
Тестовые наборы включают как степени двойки, так и произвольные размеры, что покрывает оба пути
алгоритма: только поразрядную сортировку и поразрядную сортировку с последующим слиянием Бэтчера.

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
| seq | 1 | 0.159970 | 1.00 | N/A |
| omp | 2 | 0.056191 | 2.85 | 142.34% |
| omp | 4 | 0.063745 | 2.51 | 62.73% |
| omp | 8 | 0.080437 | 1.99 | 24.86% |
| tbb | 2 | 0.123646 | 1.29 | 64.69% |
| tbb | 4 | 0.074043 | 2.16 | 54.01% |
| tbb | 8 | 0.063984 | 2.50 | 31.25% |
| stl | 2 | 0.179218 | 0.89 | 44.62% |
| stl | 4 | 0.164399 | 0.97 | 24.33% |
| stl | 8 | 0.165794 | 0.96 | 12.06% |
| all | 2 × 1 | 0.079706 | 2.01 | 100.33% |
| all | 2 × 2 | 0.099306 | 1.61 | 40.27% |
| all | 4 × 2 | 0.227350 | 0.70 | 8.79% |

OMP даёт наибольшее ускорение на двух потоках и теряет эффективность при дальнейшем росте числа потоков.
TBB масштабируется монотонно и приближается к результату OMP на восьми потоках. STL не показывает
ускорения, поскольку число рабочих потоков определяется аппаратной конфигурацией, а не заданным
параметром, и потоки пересоздаются на каждом шаге сети Бэтчера. ALL даёт ускорение в конфигурации
`2 × 1`, но проседает при росте числа процессов из-за накладных расходов на координацию MPI. Общая
причина ограниченной масштабируемости — большое число синхронизированных проходов в сети Бэтчера,
растущее логарифмически от размера массива.

## 8. Conclusions

Реализация сортирует массив `double` через поразрядную сортировку и сеть чётно-нечётного слияния Бэтчера.
Последовательная версия используется как baseline, а параллельные версии сравниваются с ней по времени
выполнения, ускорению и эффективности. Наилучшее ускорение получено у OMP на двух потоках и у ALL в
конфигурации `2 × 1`. Основное ограничение масштабируемости связано с большим числом точек синхронизации
в сети Бэтчера и накладными расходами на управление параллелизмом.

## 9. References

1. OpenMP Architecture Review Board. OpenMP Application Programming Interface.
2. oneAPI Threading Building Blocks Documentation.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: `std::thread`.

## Appendix (Optional)

Ниже приведены основные фрагменты `RunImpl()` и сети Бэтчера для всех реализаций.

### SEQ RunImpl

```cpp
bool TitaevSortirovkaBetcheraSEQ::RunImpl() {
  auto &input = GetInput();
  const size_t n = input.size();
  if (n <= 1) {
    return true;
  }

  std::vector<uint64_t> keys(n);
  ConvertToKeys(input, keys);
  RadixSort(keys);

  ConvertFromKeys(keys, GetOutput());

  if ((n & (n - 1)) == 0) {
    BatcherSort();
  }

  return true;
}
```

### OMP BatcherSort

```cpp
void TitaevSortirovkaBetcheraOMP::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
#pragma omp parallel for
      for (long long ii = 0; ii < static_cast<long long>(n); ii++) {
        const size_t i = static_cast<size_t>(ii);
        const size_t l = i ^ j;
        if (l > i) {
          const bool ascending = ((i & k) == 0);
          const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
          if (need_swap) {
            std::swap(result[i], result[l]);
          }
        }
      }
    }
  }
}
```

### TBB BatcherSort

```cpp
void TitaevSortirovkaBetcheraTBB::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
      oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, n),
                                [&](const oneapi::tbb::blocked_range<size_t> &r) {
                                  for (size_t i = r.begin(); i < r.end(); i++) {
                                    const size_t l = i ^ j;
                                    if (l > i && l < n) {
                                      const bool ascending = ((i & k) == 0);
                                      const bool need_swap =
                                          ascending ? (result[i] > result[l]) : (result[i] < result[l]);
                                      if (need_swap) {
                                        std::swap(result[i], result[l]);
                                      }
                                    }
                                  }
                                });
    }
  }
}
```

### STL BatcherSort

```cpp
void TitaevSortirovkaBetcheraSTL::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  const unsigned int num_threads = GetThreadCount();

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
      std::vector<std::thread> threads;
      threads.reserve(num_threads);
      const size_t chunk = (n + num_threads - 1) / num_threads;

      for (unsigned int t = 0; t < num_threads; t++) {
        const size_t begin = t * chunk;
        const size_t end = std::min(begin + chunk, n);
        if (begin >= end) {
          break;
        }
        threads.emplace_back([&result, n, k, j, begin, end]() {
          for (size_t i = begin; i < end; i++) {
            const size_t l = i ^ j;
            if (l > i && l < n) {
              const bool ascending = ((i & k) == 0);
              const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
              if (need_swap) {
                std::swap(result[i], result[l]);
              }
            }
          }
        });
      }
      for (auto &th : threads) {
        th.join();
      }
    }
  }
}
```

### ALL RunImpl

```cpp
bool TitaevSortirovkaBetcheraALL::RunImpl() {
  auto &input = GetInput();
  const size_t n = input.size();
  if (n <= 1) {
    return true;
  }
  std::vector<uint64_t> keys(n);
  ConvertToKeys(input, keys);
  RadixSort(keys);
  ConvertFromKeys(keys, GetOutput());
  if ((n & (n - 1)) == 0) {
    BatcherSort();
  }
  return true;
}
```
