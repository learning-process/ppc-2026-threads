# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Фокса

- Студент: Синев Артём Александрович, группа 3823Б1ПР2
- Технология: TBB
- Вариант: 2

## 1. Введение

Умножение плотных матриц является одной из фундаментальных задач параллельных вычислений. Операция активно применяется в задачах компьютерной графики, машинного обучения, научного моделирования и численных методов.

В данной работе реализован параллельный алгоритм умножения матриц с использованием библиотеки Intel oneTBB. В основе реализации лежит блочная схема алгоритма Фокса, позволяющая эффективно распределять вычисления между потоками и улучшать локальность данных.

Основной целью работы является реализация масштабируемого алгоритма умножения матриц и анализ эффективности использования TBB для многопоточных вычислений.

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

Для малых размеров матриц используется классический алгоритм умножения с тремя вложенными циклами.

```cpp
tbb::parallel_for(
    tbb::blocked_range2d<size_t>(0, n, 0, n),
    [&](const tbb::blocked_range2d<size_t> &r) {

  for (size_t i = r.rows().begin();
       i < r.rows().end();
       ++i) {

    for (size_t j = r.cols().begin();
         j < r.cols().end();
         ++j) {

      double sum = 0.0;

      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] *
               b[(k * n) + j];
      }

      c[(i * n) + j] = sum;
    }
  }
});
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

В реализации используется блочная схема алгоритма Фокса.

### Основная идея

Матрицы разбиваются на квадратные блоки размера:

```math
bs \times bs
```

Каждый поток обрабатывает отдельные блоки матриц.

### Разбиение матриц

```cpp
tbb::parallel_for(
    tbb::blocked_range2d<int>(0, q, 0, q),
    [&](const tbb::blocked_range2d<int> &r) {
```

### Выполнение шагов алгоритма Фокса

```cpp
tbb::parallel_for(
    tbb::blocked_range2d<int>(0, q, 0, q),
    [&](const tbb::blocked_range2d<int> &r) {
```

### Выбор размера блока

```cpp
size_t bs = 1;

const auto sqrt_n =
    static_cast<size_t>(
        std::sqrt(static_cast<double>(n)));

for (size_t div = sqrt_n;
     div >= 1;
     --div) {

  if (n % div == 0) {
    bs = div;
    break;
  }
}
```

### Используемые механизмы TBB

- `tbb::parallel_for`
- `tbb::blocked_range2d`
- автоматическая балансировка нагрузки
- task-based модель параллелизма

## 5. Детали реализации

### Структура проекта

- `common/include/common.hpp`
- `tbb/include/ops_tbb.hpp`
- `tbb/src/ops_tbb.cpp`
- `tests/functional/main.cpp`
- `tests/performance/main.cpp`

### Основной класс

```cpp
class SinevAMultMatrixFoxAlgorithmTBB : public BaseTask
```

### Основные методы

- `ValidationImpl()`
- `PreProcessingImpl()`
- `RunImpl()`
- `PostProcessingImpl()`

### Проверка входных данных

```cpp
bool SinevAMultMatrixFoxAlgorithmTBB::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] =
      GetInput();

  return matrix_size > 0 &&
         matrix_a.size() ==
             matrix_size * matrix_size &&
         matrix_b.size() ==
             matrix_size * matrix_size;
}
```

### Разбиение матриц на блоки

```cpp
void SinevAMultMatrixFoxAlgorithmTBB::DecomposeToBlocks(
    const std::vector<double> &src,
    std::vector<double> &dst,
    size_t n,
    size_t bs,
    int q) {

  tbb::parallel_for(
      tbb::blocked_range2d<int>(0, q, 0, q),
      [&](const tbb::blocked_range2d<int> &r) {

    for (int bi = r.rows().begin();
         bi < r.rows().end();
         ++bi) {

      for (int bj = r.cols().begin();
           bj < r.cols().end();
           ++bj) {

        const size_t block_off =
            (static_cast<size_t>((bi * q) + bj)) *
            (bs * bs);

        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {

            const size_t src_idx =
                ((static_cast<size_t>(bi) * bs + i) * n) +
                (static_cast<size_t>(bj) * bs + j);

            const size_t dst_idx =
                block_off + (i * bs) + j;

            dst[dst_idx] = src[src_idx];
          }
        }
      }
    }
  });
}
```

### Умножение блоков

```cpp
void SinevAMultMatrixFoxAlgorithmTBB::MultiplyBlocks(
    const std::vector<double> &blocks_a,
    const std::vector<double> &blocks_b,
    std::vector<double> &blocks_c,
    size_t bs,
    size_t a_off,
    size_t b_off,
    size_t c_off) {

  for (size_t ii = 0; ii < bs; ++ii) {
    for (size_t kk = 0; kk < bs; ++kk) {

      const double val =
          blocks_a[a_off + (ii * bs) + kk];

      const size_t b_base =
          b_off + (kk * bs);

      const size_t c_base =
          c_off + (ii * bs);

      for (size_t jj = 0; jj < bs; ++jj) {

        blocks_c[c_base + jj] +=
            val *
            blocks_b[b_base + jj];
      }
    }
  }
}
```

### Выполнение шага алгоритма Фокса

```cpp
void SinevAMultMatrixFoxAlgorithmTBB::FoxStep(
    const std::vector<double> &blocks_a,
    const std::vector<double> &blocks_b,
    std::vector<double> &blocks_c,
    size_t bs,
    int q,
    int step) {

  const size_t block_size = bs * bs;

  tbb::parallel_for(
      tbb::blocked_range2d<int>(0, q, 0, q),
      [&](const tbb::blocked_range2d<int> &r) {

    for (int i = r.rows().begin();
         i < r.rows().end();
         ++i) {

      for (int j = r.cols().begin();
           j < r.cols().end();
           ++j) {

        const int k = (i + step) % q;

        const size_t a_off =
            (static_cast<size_t>((i * q) + k)) *
            block_size;

        const size_t b_off =
            (static_cast<size_t>((k * q) + j)) *
            block_size;

        const size_t c_off =
            (static_cast<size_t>((i * q) + j)) *
            block_size;

        MultiplyBlocks(
            blocks_a,
            blocks_b,
            blocks_c,
            bs,
            a_off,
            b_off,
            c_off);
      }
    }
  });
}
```

### Особенности реализации

- Используется task-based модель TBB
- Применяется двумерное разбиение диапазонов
- Используется автоматическая балансировка нагрузки
- Улучшается локальность данных
- Для малых матриц используется простое умножение

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
./build/bin/ppc_func_tests --gtest_filter="*Sinev*tbb*" -v
```

Запуск performance тестов:

```bash
./build/bin/ppc_perf_tests --gtest_filter="*Sinev*tbb*" -v
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

Ниже приведено сравнение реализаций алгоритма Фокса для последовательной версии, OpenMP и Intel TBB.

| Режим              | Потоки | Время, с   | Ускорение | Эффективность |
|--------------------|---------|-------------|------------|----------------|
| seq                | 1       | 0.0261574700 | 1.00       | 100%           |
| omp (pipeline)     | 16      | 0.0024294812 | 10.77      | 67.3%          |
| omp (task_run)     | 16      | 0.0042057726 | 6.22       | 38.9%          |
| tbb (pipeline)     | 16      | 0.0166006096 | 1.58       | 9.9%           |
| tbb (task_run)     | 16      | 0.0166640150 | 1.57       | 9.8%           |

Ускорение:

```math
Speedup = \frac{T_{seq}}{T_{omp}}
```

Эффективность:

```math
Efficiency = \frac{Speedup}{Threads} \times 100\%
```

#### Анализ результатов

Реализация Intel TBB показала меньший прирост производительности. Несмотря на использование 16 потоков, ускорение составило около 1.6 раза. Основными причинами являются:

- высокие накладные расходы TBB на создание и управление задачами;
- малый размер тестовых матриц;
- частые обращения к памяти при работе с блоками;
- недостаточная вычислительная нагрузка для полного насыщения всех потоков.

Таким образом:

- OpenMP оказался наиболее эффективной технологией для данной задачи;
- TBB обеспечивает корректное распараллеливание, однако уступает OpenMP по скорости;
- последовательная реализация значительно проигрывает параллельным версиям при увеличении размера матриц.

При увеличении размерности матриц ожидается рост эффективности TBB за счёт лучшей загрузки потоков и уменьшения относительных накладных расходов.

### Узкие места алгоритма

Основные ограничения производительности:

1. Накладные расходы task scheduler
2. Конкуренция за кэш-память
3. Большое количество операций записи
4. Ограничение пропускной способности памяти

### Масштабируемость

TBB обеспечивает:

- автоматическое распределение задач
- динамическую балансировку нагрузки
- эффективную работу с многопоточностью
- снижение накладных расходов ручного управления потоками

## 8. Выводы

В ходе работы была реализована TBB версия алгоритма Фокса для умножения плотных матриц.

### Полученные результаты

- Реализован параллельный алгоритм на основе Intel oneTBB
- Успешно пройдены все функциональные тесты
- Получено ускорение относительно последовательной версии
- Реализована блочная схема вычислений

### Особенности реализации

- Использование task-based параллелизма
- Автоматическая балансировка нагрузки
- Улучшенная локальность данных
- Эффективное распределение блоков

### Ограничения

- Производительность зависит от размера блоков
- Возможны накладные расходы scheduler
- Ограничение масштабирования памятью

### Перспективы развития

- Автоматический подбор оптимального размера блока
- Использование SIMD-инструкций
- NUMA-оптимизации
- Гибридные схемы MPI + TBB

## 9. Источники

1. Лекции по параллельному программированию Сысоева А. В.
2. Документация Intel oneTBB: <https://oneapi-src.github.io/oneTBB/>
3. Материалы курса: <https://github.com/learning-process/ppc-2026-threads>
4. Fox G. C. Matrix Algorithms on Parallel Hardware.
5. oneAPI Threading Building Blocks Developer Guide

## 10. Приложение

```cpp
#include "sinev_a_mult_matrix_fox_algorithm/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmTBB::SinevAMultMatrixFoxAlgorithmTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmTBB::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmTBB::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

// Добавляем static к определениям
void SinevAMultMatrixFoxAlgorithmTBB::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
  tbb::parallel_for(tbb::blocked_range2d<size_t>(0, n, 0, n), [&](const tbb::blocked_range2d<size_t> &r) {
    for (size_t i = r.rows().begin(); i < r.rows().end(); ++i) {
      for (size_t j = r.cols().begin(); j < r.cols().end(); ++j) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        c[(i * n) + j] = sum;
      }
    }
  });
}

void SinevAMultMatrixFoxAlgorithmTBB::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, q, 0, q), [&](const tbb::blocked_range2d<int> &r) {
    for (int bi = r.rows().begin(); bi < r.rows().end(); ++bi) {
      for (int bj = r.cols().begin(); bj < r.cols().end(); ++bj) {
        const size_t block_off = (static_cast<size_t>((bi * q) + bj)) * (bs * bs);
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
            const size_t dst_idx = block_off + (i * bs) + j;
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    }
  });
}

void SinevAMultMatrixFoxAlgorithmTBB::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, q, 0, q), [&](const tbb::blocked_range2d<int> &r) {
    for (int bi = r.rows().begin(); bi < r.rows().end(); ++bi) {
      for (int bj = r.cols().begin(); bj < r.cols().end(); ++bj) {
        const size_t block_off = (static_cast<size_t>((bi * q) + bj)) * (bs * bs);
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            const size_t src_idx = block_off + (i * bs) + j;
            const size_t dst_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    }
  });
}

void SinevAMultMatrixFoxAlgorithmTBB::MultiplyBlocks(const std::vector<double> &blocks_a,
                                                     const std::vector<double> &blocks_b, std::vector<double> &blocks_c,
                                                     size_t bs, size_t a_off, size_t b_off, size_t c_off) {
  for (size_t ii = 0; ii < bs; ++ii) {
    for (size_t kk = 0; kk < bs; ++kk) {
      const double val = blocks_a[a_off + (ii * bs) + kk];
      const size_t b_base = b_off + (kk * bs);
      const size_t c_base = c_off + (ii * bs);
      for (size_t jj = 0; jj < bs; ++jj) {
        blocks_c[c_base + jj] += val * blocks_b[b_base + jj];
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmTBB::FoxStep(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                                              std::vector<double> &blocks_c, size_t bs, int q, int step) {
  const size_t block_size = bs * bs;

  tbb::parallel_for(tbb::blocked_range2d<int>(0, q, 0, q), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
        const int k = (i + step) % q;

        const size_t a_off = (static_cast<size_t>((i * q) + k)) * block_size;
        const size_t b_off = (static_cast<size_t>((k * q) + j)) * block_size;
        const size_t c_off = (static_cast<size_t>((i * q) + j)) * block_size;

        MultiplyBlocks(blocks_a, blocks_b, blocks_c, bs, a_off, b_off, c_off);
      }
    }
  });
}

bool SinevAMultMatrixFoxAlgorithmTBB::RunImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);
  auto &c = GetOutput();

  // Для маленьких матриц используем простое умножение
  if (n <= 8) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  size_t bs = ChooseBlockSize(n);
  const int actual_q = static_cast<int>(n / bs);

  const auto total_blocks = static_cast<size_t>(actual_q) * static_cast<size_t>(actual_q);
  const auto block_elements = bs * bs;

  std::vector<double> blocks_a(total_blocks * block_elements);
  std::vector<double> blocks_b(total_blocks * block_elements);
  std::vector<double> blocks_c(total_blocks * block_elements, 0.0);

  DecomposeToBlocks(a, blocks_a, n, bs, actual_q);
  DecomposeToBlocks(b, blocks_b, n, bs, actual_q);

  for (int step = 0; step < actual_q; ++step) {
    FoxStep(blocks_a, blocks_b, blocks_c, bs, actual_q, step);
  }

  AssembleFromBlocks(blocks_c, c, n, bs, actual_q);

  return true;
}

size_t SinevAMultMatrixFoxAlgorithmTBB::ChooseBlockSize(size_t n) {
  size_t bs = 1;
  const auto sqrt_n = static_cast<size_t>(std::sqrt(static_cast<double>(n)));
  for (size_t div = sqrt_n; div >= 1; --div) {
    if (n % div == 0) {
      bs = div;
      break;
    }
  }
  return bs;
}

bool SinevAMultMatrixFoxAlgorithmTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm

```
