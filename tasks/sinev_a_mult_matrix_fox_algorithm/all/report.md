# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Фокса

- Студент: Синев Артём Александрович, группа 3823Б1ПР2
- Технология: MPI+OpenMP (ALL)
- Вариант: 2

## 1. Введение

Умножение плотных матриц является одной из фундаментальных задач параллельных вычислений. Операция активно применяется в задачах компьютерной графики, машинного обучения, научного моделирования и численных методов.

В данной работе реализован гибридный параллельный алгоритм умножения матриц с использованием MPI (Message Passing Interface) для межпроцессного взаимодействия и OpenMP для внутрипроцессного распараллеливания. В основе реализации лежит блочная схема алгоритма Фокса, позволяющая эффективно распределять вычисления между процессами и потоками.

Основной целью работы является реализация масштабируемого гибридного алгоритма умножения матриц с использованием MPI+OpenMP и анализ эффективности по сравнению с предыдущими реализациями (STL, OpenMP, Intel TBB).

## 2. Постановка задачи

Даны две квадратные матрицы `A` и `B` размера:

```math
N \times N
```

Необходимо вычислить произведение:

```math
C = A \times B
```

где:

```math
C_{ij} =
\sum_{k=0}^{N-1}
A_{ik} \cdot B_{kj}
```

### Входные данные

- `matrix_size` — размер матрицы
- `matrix_a` — первая матрица
- `matrix_b` — вторая матрица

### Выходные данные

- `output` — результирующая матрица

### Ограничения

- Размер матриц должен быть больше нуля
- Матрицы должны быть квадратными
- Размер входных массивов должен соответствовать `N × N`
- Тип данных элементов — `double`

## 3. Описание базового алгоритма

Для случая fallback (когда количество процессов не является квадратом или размер матрицы не кратен количеству процессов) используется классический алгоритм умножения с тремя вложенными циклами, распараллеленный через OpenMP.

```cpp
#pragma omp parallel for default(none) shared(n, a, b, c) collapse(2)
for (size_t i = 0; i < n; ++i) {
  for (size_t j = 0; j < n; ++j) {
    double sum = 0.0;
    for (size_t k = 0; k < n; ++k) {
      sum += a[(i * n) + k] * b[(k * n) + j];
    }
    c[(i * n) + j] = sum;
  }
}
```

### Сложность алгоритма

#### Временная сложность

```math
O(N^3)
```

#### Пространственная сложность

```math
O(N^2)
```

## 4. Схема распараллеливания

В реализации используется гибридная схема: MPI для распределения блоков между процессами и OpenMP для параллельной обработки внутри каждого процесса.

### Основная идея

Матрицы разбиваются на квадратные блоки размера:

```math
bs \times bs
```

Количество процессов `world_size` должно быть точным квадратом: `q = sqrt(world_size)`. Процессы организуются в двумерную решетку `q × q`. Каждый процесс получает один блок матрицы A и один блок матрицы B.

### MPI-коммуникации

- `MPI_Scatter` — распределение блоков от root-процесса всем процессам
- `MPI_Bcast` — широковещательная рассылка блоков A внутри строки процессов
- `MPI_Sendrecv_replace` — циклический сдвиг блоков B между процессами
- `MPI_Gather` — сбор результатов от всех процессов
- `MPI_Comm_split` — создание коммуникаторов для строк решетки

### Внутрипроцессное распараллеливание

```cpp
#pragma omp parallel for default(none) shared(local_a, local_b, local_c, bs) collapse(2)
for (size_t i = 0; i < bs; ++i) {
  for (size_t j = 0; j < bs; ++j) {
    double sum = 0.0;
    for (size_t k = 0; k < bs; ++k) {
      sum += local_a[(i * bs) + k] * local_b[(k * bs) + j];
    }
    local_c[(i * bs) + j] += sum;
  }
}
```

### Используемые механизмы

- `MPI_Init` / `MPI_Finalize` — инициализация MPI
- `MPI_Comm_rank` / `MPI_Comm_size` — идентификация процессов
- `MPI_Scatter` / `MPI_Gather` — распределение и сбор данных
- `MPI_Bcast` — широковещательная рассылка
- `MPI_Sendrecv_replace` — коммуникация с заменой буфера
- `MPI_Comm_split` — создание подкоммуникаторов
- `#pragma omp parallel for collapse(2)` — параллельные циклы OpenMP

## 5. Детали реализации

### Структура проекта

- `common/include/common.hpp`
- `all/include/ops_all.hpp`
- `all/src/ops_all.cpp`
- `tests/functional/main.cpp`
- `tests/performance/main.cpp`

### Основной класс

```cpp
class SinevAMultMatrixFoxAlgorithmALL : public BaseTask
```

### Основные методы

- `ValidationImpl()`
- `PreProcessingImpl()`
- `RunImpl()`
- `PostProcessingImpl()`

### Проверка входных данных

```cpp
bool SinevAMultMatrixFoxAlgorithmALL::ValidationImpl() {
  const auto &[n, a, b] = GetInput();
  return (n > 0U) && (a.size() == (n * n)) && (b.size() == (n * n));
}
```

### Разбиение матриц на блоки

```cpp
void SinevAMultMatrixFoxAlgorithmALL::DecomposeToBlocks(
    const std::vector<double> &src,
    std::vector<double> &dst,
    size_t n, size_t bs, int q) {
  
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_offset = static_cast<size_t>((bi * q) + bj) * (bs * bs);
      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + 
                                 (static_cast<size_t>(bj) * bs + j);
          const size_t dst_idx = block_offset + (i * bs) + j;
          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}
```

### Умножение локальных блоков

```cpp
void SinevAMultMatrixFoxAlgorithmALL::LocalMatrixMultiply(
    const std::vector<double> &local_a,
    const std::vector<double> &local_b,
    std::vector<double> &local_c, size_t bs) {
  
#pragma omp parallel for default(none) shared(local_a, local_b, local_c, bs) collapse(2)
  for (size_t i = 0; i < bs; ++i) {
    for (size_t j = 0; j < bs; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < bs; ++k) {
        sum += local_a[(i * bs) + k] * local_b[(k * bs) + j];
      }
      local_c[(i * bs) + j] += sum;
    }
  }
}
```

### Алгоритм Фокса с MPI

```cpp
void SinevAMultMatrixFoxAlgorithmALL::RunFoxStages(
    int q, int row, int col, size_t bs, size_t block_size,
    MPI_Comm row_comm, std::vector<double> &local_a,
    std::vector<double> &local_b, std::vector<double> &local_c) {
  
  std::vector<double> temp_a(block_size);

  for (int step = 0; step < q; ++step) {
    const int root = (row + step) % q;

    if (col == root) {
      temp_a = local_a;
    }

    MPI_Bcast(temp_a.data(), static_cast<int>(block_size), MPI_DOUBLE, root, row_comm);

    LocalMatrixMultiply(temp_a, local_b, local_c, bs);

    const int send_to = (((row - 1 + q) % q) * q) + col;
    const int recv_from = (((row + 1) % q) * q) + col;

    MPI_Sendrecv_replace(local_b.data(), static_cast<int>(block_size), MPI_DOUBLE,
                         send_to, 0, recv_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}
```

### Основной вычислительный метод

```cpp
bool SinevAMultMatrixFoxAlgorithmALL::RunImpl() {
  int rank = 0, world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto &[n, a, b] = GetInput();
  auto &c = GetOutput();

  const int q = static_cast<int>(std::sqrt(world_size));

  if (NeedFallback(n, q, world_size)) {
    ExecuteFallback(rank, n, a, b, c);
    return true;
  }

  const size_t bs = n / static_cast<size_t>(q);
  const size_t block_size = bs * bs;
  const int row = rank / q;
  const int col = rank % q;

  std::vector<double> local_a(block_size);
  std::vector<double> local_b(block_size);
  std::vector<double> local_c(block_size, 0.0);

  std::vector<double> blocks_a, blocks_b;

  if (rank == 0) {
    blocks_a.resize(static_cast<size_t>(world_size) * block_size);
    blocks_b.resize(static_cast<size_t>(world_size) * block_size);
    DecomposeToBlocks(a, blocks_a, n, bs, q);
    DecomposeToBlocks(b, blocks_b, n, bs, q);
  }

  ScatterBlocks(rank, blocks_a, blocks_b, local_a, local_b, block_size);

  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, row, col, &row_comm);

  RunFoxStages(q, row, col, bs, block_size, row_comm, local_a, local_b, local_c);

  GatherResult(rank, world_size, n, bs, block_size, q, local_c, c);

  if (row_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&row_comm);
  }

  return true;
}
```

### Особенности реализации

- Гибридная модель MPI+OpenMP для максимальной производительности
- Автоматический fallback при некорректном количестве процессов
- Использование `MPI_Sendrecv_replace` для эффективной коммуникации
- OpenMP collapse(2) для улучшения локальности данных
- Создание строковых коммуникаторов через `MPI_Comm_split`

## 6. Экспериментальное окружение

### 6.1 Аппаратное обеспечение / ОС

- **Процессор:** Intel Core i7-13700HX
- **Количество ядер:** 16
- **ОЗУ:** 8 ГБ
- **ОС:** Kubuntu 24.04

### 6.2 Программное окружение

- **Компилятор:** g++ 13.3.0
- **Стандарт C++:** C++20
- **Библиотека:** Intel oneTBB
- **Тип сборки:** Release
- **Система сборки:** CMake

### 6.3 Тестовое окружение

Запуск функциональных тестов:

```bash
mpirun -n 4 ./build/bin/ppc_func_tests --gtest_filter="*Sinev*all*" -v
```

Запуск performance тестов:

```bash
mpirun -n <N> ./build/bin/ppc_perf_tests --gtest_filter="*Sinev*all*" -v
```

## 7. Результаты

### 7.1 Корректность работы

Были выполнены тесты для матриц размеров:

```text
1×1
2×2
4×4
6×6
8×8
9×9
10×10
16×16
18×18
25×25
50×50
75×75
100×100
```

Все тесты завершились успешно.

```text
[==========] Running 13 tests from 1 test suite.
[  PASSED  ] 13 tests.
```

### Проверка корректности

Проверка выполнялась сравнением с эталонной последовательной реализацией.

Используемая точность сравнения:

```cpp
const double epsilon = 1e-10;
```

### 7.2 Производительность

Ниже приведено сравнение реализаций алгоритма Фокса для последовательной версии, OpenMP, Intel TBB и STL.

| Режим        | Процессы      | Потоки | Время, с   | Ускорение | Эффективность |
|--------------|--------------|---------|-------------|------------|----------------|
| seq           |      1       | 1       | 0.0261574700 | 1.00       | 100%           |
| omp (pipeline)|     1        | 16      | 0.0024294812 | 10.77      | 67.3%          |
| tbb (pipeline)|      1       | 16      | 0.0166006096 | 1.58       | 9.9%           |
| stl (pipeline)|       1     | 16      | 0.0124331014 | 2.10       | 13.1%          |
| MPI+OpenMP|       4     | 4      | 0.004737| 5.52       | 34.5%           |
| MPI+OpenMP|        5    | -      | 0.007615 | 3.44      | 17.2%          |
| MPI+OpenMP|        8     | 2      | 0.009797 | 2.67       | 8.3%          |
| MPI+OpenMP|        12      | -      | 0.018164 | 1.44     | 3.0%          |

Ускорение:

```math
Speedup = \frac{T_{seq}}{T_{omp}}
```

Эффективность:

```math
Efficiency = \frac{Speedup}{Threads} \times 100\%
```

#### Анализ результатов

MPI+OpenMP (гибридная реализация):

- Максимальное ускорение достигнуто при 4 процессах (по 4 потока в каждом) — 5.52x
- При 2 процессах (по 8 потоков) ускорение составляет 2.67x
- При 3 процессах (fallback) — 3.44x
- При 4 процессах эффективность составила 34.5%
- При увеличении числа процессов сверх оптимального значения производительность падает из-за накладных расходов на коммуникации

### Основные причины различий

1. OpenMP остается лидером для одноузловых вычислений благодаря минимальным накладным расходам и эффективному использованию кэш-памяти

2. MPI+OpenMP показывает лучший баланс:

- Ускорение 5.52x при 4 процессах (16 потоков всего)
- Позволяет масштабироваться на многоузловые системы
- Накладные расходы на коммуникации ограничивают масштабируемость

3. STL и TBB уступают специализированным решениям из-за:

- Отсутствия оптимизированных коллективных операций
- Ручного управления распределением данных

### Узкие места алгоритма

- Накладные расходы на MPI-коммуникации при каждом шаге алгоритма
- Требование к точному квадрату количества процессов
- Необходимость синхронизации между шагами алгоритма Фокса
- Дополнительные накладные расходы на разбиение и сборку блоков

### Масштабируемость

MPI+OpenMP обеспечивает:

- Возможность распределения вычислений на несколько узлов
- Балансировку нагрузки через блочное разбиение
- Эффективное использование многоядерных процессоров внутри узла
- Масштабирование до 4-8 процессов на данной конфигурации

## 8. Выводы

В ходе работы была реализована гибридная MPI+OpenMP версия алгоритма Фокса для умножения плотных матриц.

### Полученные результаты

- Реализован гибридный параллельный алгоритм на основе MPI+OpenMP
- Успешно пройдены все функциональные тесты (13 тестов)
- Получено ускорение до 5.52x относительно последовательной версии
- Реализован fallback-механизм для некорректных конфигураций
- Достигнута эффективность 34.5% при 4 процессах

### Особенности реализации

- Гибридная модель позволяет использовать преимущества обоих подходов
- MPI обеспечивает распределение данных и коммуникации
- OpenMP обеспечивает эффективные вычисления внутри процесса
- Автоматический fallback для некорректных конфигураций

### Ограничения

- Количество процессов должно быть точным квадратом
- Размер матрицы должен делиться на количество процессов в строке/столбце
- Накладные расходы на MPI-коммуникации ограничивают масштабируемость
- Требуется сбалансированное распределение блоков

### Перспективы развития

- Использование неблокирующих MPI-коммуникаций
- Реализация динамической балансировки нагрузки
- Оптимизация размера блока под кэш-иерархию
- Поддержка неквадратных решеток процессов
- Использование MPI-внутренних оптимизаций коллективных операций

## 9. Источники

1. Лекции по параллельному программированию Сысоева А. В.
2. Материалы курса: <https://github.com/learning-process/ppc-2026-threads>
3. Fox G. C. Matrix Algorithms on Parallel Hardware

## 10. Приложение

```cpp
#include "sinev_a_mult_matrix_fox_algorithm/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmALL::SinevAMultMatrixFoxAlgorithmALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmALL::ValidationImpl() {
  const auto &[n, a, b] = GetInput();

  return (n > 0U) && (a.size() == (n * n)) && (b.size() == (n * n));
}

bool SinevAMultMatrixFoxAlgorithmALL::PreProcessingImpl() {
  const auto &[n, a, b] = GetInput();

  GetOutput().resize(n * n, 0.0);

  return true;
}

void SinevAMultMatrixFoxAlgorithmALL::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
#pragma omp parallel for default(none) shared(n, a, b, c) collapse(2)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;

      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }

      c[(i * n) + j] = sum;
    }
  }
}

void SinevAMultMatrixFoxAlgorithmALL::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_offset = static_cast<size_t>((bi * q) + bj) * (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);

          const size_t dst_idx = block_offset + (i * bs) + j;

          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmALL::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_offset = static_cast<size_t>((bi * q) + bj) * (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = block_offset + (i * bs) + j;

          const size_t dst_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);

          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmALL::LocalMatrixMultiply(const std::vector<double> &local_a,
                                                          const std::vector<double> &local_b,
                                                          std::vector<double> &local_c, size_t bs) {
#pragma omp parallel for default(none) shared(local_a, local_b, local_c, bs) collapse(2)
  for (size_t i = 0; i < bs; ++i) {
    for (size_t j = 0; j < bs; ++j) {
      double sum = 0.0;

      for (size_t k = 0; k < bs; ++k) {
        sum += local_a[(i * bs) + k] * local_b[(k * bs) + j];
      }

      local_c[(i * bs) + j] += sum;
    }
  }
}

bool SinevAMultMatrixFoxAlgorithmALL::NeedFallback(size_t n, int q, int world_size) {
  return ((q * q) != world_size) || ((n % static_cast<size_t>(q)) != 0U);
}

void SinevAMultMatrixFoxAlgorithmALL::ExecuteFallback(int rank, size_t n, const std::vector<double> &a,
                                                      const std::vector<double> &b, std::vector<double> &c) {
  if (rank == 0) {
    SimpleMultiply(n, a, b, c);
  }

  MPI_Bcast(c.data(), static_cast<int>(n * n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void SinevAMultMatrixFoxAlgorithmALL::ScatterBlocks(int rank, const std::vector<double> &blocks_a,
                                                    const std::vector<double> &blocks_b, std::vector<double> &local_a,
                                                    std::vector<double> &local_b, size_t block_size) {
  const double *send_a = (rank == 0) ? blocks_a.data() : nullptr;

  const double *send_b = (rank == 0) ? blocks_b.data() : nullptr;

  MPI_Scatter(send_a, static_cast<int>(block_size), MPI_DOUBLE, local_a.data(), static_cast<int>(block_size),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Scatter(send_b, static_cast<int>(block_size), MPI_DOUBLE, local_b.data(), static_cast<int>(block_size),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void SinevAMultMatrixFoxAlgorithmALL::RunFoxStages(int q, int row, int col, size_t bs, size_t block_size,
                                                   MPI_Comm row_comm, std::vector<double> &local_a,
                                                   std::vector<double> &local_b, std::vector<double> &local_c) {
  std::vector<double> temp_a(block_size);

  for (int step = 0; step < q; ++step) {
    const int root = (row + step) % q;

    if (col == root) {
      temp_a = local_a;
    }

    MPI_Bcast(temp_a.data(), static_cast<int>(block_size), MPI_DOUBLE, root, row_comm);

    LocalMatrixMultiply(temp_a, local_b, local_c, bs);

    const int send_to = (((row - 1 + q) % q) * q) + col;

    const int recv_from = (((row + 1) % q) * q) + col;

    MPI_Sendrecv_replace(local_b.data(), static_cast<int>(block_size), MPI_DOUBLE, send_to, 0, recv_from, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void SinevAMultMatrixFoxAlgorithmALL::GatherResult(int rank, int world_size, size_t n, size_t bs, size_t block_size,
                                                   int q, const std::vector<double> &local_c, std::vector<double> &c) {
  std::vector<double> blocks_c;

  if (rank == 0) {
    blocks_c.resize(static_cast<size_t>(world_size) * block_size);
  }

  double *recv_buffer = (rank == 0) ? blocks_c.data() : nullptr;

  MPI_Gather(local_c.data(), static_cast<int>(block_size), MPI_DOUBLE, recv_buffer, static_cast<int>(block_size),
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    AssembleFromBlocks(blocks_c, c, n, bs, q);
  }

  MPI_Bcast(c.data(), static_cast<int>(n * n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

bool SinevAMultMatrixFoxAlgorithmALL::RunImpl() {
  int rank = 0;
  int world_size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto &[n, a, b] = GetInput();

  auto &c = GetOutput();

  const int q = static_cast<int>(std::sqrt(world_size));

  if (NeedFallback(n, q, world_size)) {
    ExecuteFallback(rank, n, a, b, c);

    return true;
  }

  const size_t bs = n / static_cast<size_t>(q);

  const size_t block_size = bs * bs;

  const int row = rank / q;
  const int col = rank % q;

  std::vector<double> local_a(block_size);
  std::vector<double> local_b(block_size);
  std::vector<double> local_c(block_size, 0.0);

  std::vector<double> blocks_a;
  std::vector<double> blocks_b;

  if (rank == 0) {
    blocks_a.resize(static_cast<size_t>(world_size) * block_size);

    blocks_b.resize(static_cast<size_t>(world_size) * block_size);

    DecomposeToBlocks(a, blocks_a, n, bs, q);

    DecomposeToBlocks(b, blocks_b, n, bs, q);
  }

  ScatterBlocks(rank, blocks_a, blocks_b, local_a, local_b, block_size);

  MPI_Comm row_comm = MPI_COMM_NULL;

  const int color = row;
  const int key = col;

  MPI_Comm_split(MPI_COMM_WORLD, color, key, &row_comm);

  RunFoxStages(q, row, col, bs, block_size, row_comm, local_a, local_b, local_c);

  GatherResult(rank, world_size, n, bs, block_size, q, local_c, c);

  if (row_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&row_comm);
  }

  return true;
}

bool SinevAMultMatrixFoxAlgorithmALL::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm

```
