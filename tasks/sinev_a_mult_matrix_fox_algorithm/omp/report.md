# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Фокса.

- Студент: Синев Артём Александрович, группа 3823Б1ПР2
- Технология: OMP
- Вариант: 2

## 1. Введение

Умножение плотных матриц является одной из ключевых задач высокопроизводительных вычислений. Данная операция активно используется в линейной алгебре, машинном обучении, моделировании физических процессов и научных вычислениях.

В данной работе реализован параллельный алгоритм умножения квадратных матриц с использованием OpenMP. В основе реализации лежит блочная схема алгоритма Фокса, позволяющая эффективно распределять вычисления между потоками и улучшать локальность данных.

Цель работы — реализовать многопоточную версию алгоритма умножения матриц, проверить корректность вычислений и оценить производительность OpenMP реализации.

## 2. Постановка задачи

Даны две квадратные плотные матрицы `A` и `B` размера `N × N`, содержащие элементы типа `double`.

Требуется вычислить произведение матриц:

```math
C = A \times B
```

где:

```math
C_{ij} = \sum_{k=0}^{N-1} A_{ik} \cdot B_{kj}
```

### Входные данные

- `matrix_size` — размер матрицы `N`
- `matrix_a` — первая матрица размера `N × N`
- `matrix_b` — вторая матрица размера `N × N`

### Выходные данные

- `output` — результирующая матрица размера `N × N`

### Ограничения

- Размер матрицы должен быть больше нуля
- Матрицы должны быть квадратными
- Размер входных массивов должен соответствовать `N × N`
- Тип элементов — `double`

## 3. Описание базового алгоритма

В качестве основы используется классический алгоритм умножения матриц с тремя вложенными циклами.

Для малых матриц (`N <= 8`) используется прямое умножение:

```cpp
#pragma omp parallel for default(none) shared(n, a, b, c) collapse(2)
for (size_t i = 0; i < n; ++i) {
  for (size_t j = 0; j < n; ++j) {
    double sum = 0.0;

    for (size_t k = 0; k < n; ++k) {
      sum += a[(i * n) + k] *
             b[(k * n) + j];
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

В реализации используется блочная схема алгоритма Фокса.

### Основная идея

Матрицы разбиваются на блоки размера:

```math
bs \times bs
```

После чего:
- блоки распределяются между потоками
- каждый поток вычисляет произведение соответствующих блоков
- результат записывается в локальный блок результирующей матрицы

### Выбор размера блока

Размер блока вычисляется автоматически:

```cpp
size_t bs = 1;

auto sqrt_n =
    static_cast<size_t>(std::sqrt(static_cast<double>(n)));

for (size_t div = sqrt_n; div >= 1; --div) {
  if (n % div == 0) {
    bs = div;
    break;
  }
}
```

### Разбиение матриц на блоки

```cpp
#pragma omp parallel for default(none) \
shared(src, dst, n, bs, q) collapse(2)
for (int bi = 0; bi < q; ++bi) {
  for (int bj = 0; bj < q; ++bj) {
```

### Параллельное выполнение шагов Фокса

```cpp
#pragma omp parallel for default(none) \
shared(blocks_a, blocks_b, blocks_c,
       bs, q, step, block_size_bytes) collapse(2)
for (int i = 0; i < q; ++i) {
  for (int j = 0; j < q; ++j) {
```

### Используемые механизмы OpenMP

- `parallel for`
- `collapse(2)`
- разделяемая память (`shared`)
- автоматическое распределение итераций

## 5. Детали реализации

### Структура проекта

- `common/include/common.hpp` — общие типы данных
- `omp/include/ops_omp.hpp` — объявление OpenMP версии
- `omp/src/ops_omp.cpp` — реализация алгоритма
- `tests/functional/main.cpp` — функциональные тесты
- `tests/performance/main.cpp` — тесты производительности

### Основной класс

```cpp
class SinevAMultMatrixFoxAlgorithmOMP : public BaseTask
```

### Основные методы

- `ValidationImpl()`
- `PreProcessingImpl()`
- `RunImpl()`
- `PostProcessingImpl()`

### Проверка входных данных

```cpp
bool SinevAMultMatrixFoxAlgorithmOMP::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 &&
         matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}
```

### Разбиение матриц на блоки

```cpp
void SinevAMultMatrixFoxAlgorithmOMP::DecomposeToBlocks(
    const std::vector<double> &src,
    std::vector<double> &dst,
    size_t n,
    size_t bs,
    int q) {

#pragma omp parallel for default(none) \
shared(src, dst, n, bs, q) collapse(2)

  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
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
}
```

### Шаг алгоритма Фокса

```cpp
void SinevAMultMatrixFoxAlgorithmOMP::FoxStep(
    const std::vector<double> &blocks_a,
    const std::vector<double> &blocks_b,
    std::vector<double> &blocks_c,
    size_t bs,
    int q,
    int step) {

  const size_t block_size_bytes = bs * bs;

#pragma omp parallel for default(none) \
shared(blocks_a, blocks_b, blocks_c,
       bs, q, step, block_size_bytes) collapse(2)

  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) {

      const int k = (i + step) % q;

      const size_t a_off =
          (static_cast<size_t>((i * q) + k)) *
          block_size_bytes;

      const size_t b_off =
          (static_cast<size_t>((k * q) + j)) *
          block_size_bytes;

      const size_t c_off =
          (static_cast<size_t>((i * q) + j)) *
          block_size_bytes;

      for (size_t ii = 0; ii < bs; ++ii) {
        for (size_t kk = 0; kk < bs; ++kk) {

          const double val =
              blocks_a[a_off + (ii * bs) + kk];

          for (size_t jj = 0; jj < bs; ++jj) {

            blocks_c[c_off + (ii * bs) + jj] +=
                val *
                blocks_b[b_off + (kk * bs) + jj];
          }
        }
      }
    }
  }
}
```

### Особенности реализации

- Используется блочное представление матриц
- Применяется OpenMP распараллеливание
- Используется `collapse(2)` для лучшей балансировки
- Для малых матриц используется обычное умножение
- Алгоритм оптимизирован под кэш-локальность

## 6. Экспериментальное окружение

### 6.1 Аппаратное обеспечение / ОС

- **Процессор:** Intel Core i7-13700HX
- **Количество ядер:** 16
- **ОЗУ:** 8 ГБ
- **ОС:** Kubuntu 24.04

### 6.2 Программное окружение

- **Компилятор:** g++ 13.3.0
- **Стандарт C++:** C++20
- **OpenMP:** libgomp
- **Тип сборки:** Release
- **Система сборки:** CMake

### 6.3 Тестовое окружение

Запуск функциональных тестов:

```bash
./build/bin/ppc_func_tests --gtest_filter="*Sinev*omp*" -v
```

Запуск performance тестов:

```bash
./build/bin/ppc_perf_tests --gtest_filter="*Sinev*omp*" -v
```

## 7. Результаты

### 7.1 Корректность работы

Были выполнены функциональные тесты для различных размеров матриц:

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

Все тесты были успешно пройдены.

Результат запуска функциональных тестов:

```text
[==========] Running 13 tests from 1 test suite.
[  PASSED  ] 13 tests.
```

### Проверка результатов

Для проверки корректности использовалась эталонная последовательная реализация.

Сравнение результатов выполнялось с точностью:

```cpp
const double epsilon = 1e-10;
```

### 7.2 Производительность

Были выполнены performance-тесты OpenMP версии алгоритма.

| Режим | Количество потоков | Время, с | Ускорение | Эффективность |
|------|-------------------|----------|------------|----------------|
| seq | 1 | 0.0261574700 | 1.00 | 100% |
| omp (pipeline) | 16 | 0.0024294812 | 10.77 | 67.3% |
| omp (task_run) | 16 | 0.0042057726 | 6.22 | 38.9% |

### Формулы вычисления

Ускорение:

```math
Speedup = \frac{T_{seq}}{T_{omp}}
```

Эффективность:

```math
Efficiency = \frac{Speedup}{Threads} \times 100\%
```

### Анализ производительности

- OpenMP версия демонстрирует значительное ускорение относительно последовательной реализации
- Наилучший результат показал режим `pipeline`
- Использование блочной схемы улучшает локальность данных
- Распараллеливание циклов позволяет эффективно задействовать ядра процессора

### Узкие места алгоритма

Основными ограничениями являются:

1. Накладные расходы OpenMP
2. Конкуренция за кэш-память
3. Ограниченная пропускная способность памяти
4. Неравномерность нагрузки для некоторых размеров блоков

### Масштабируемость

Алгоритм хорошо масштабируется благодаря:
- независимости вычислений блоков
- эффективному распределению работы
- использованию `collapse(2)`
- минимальному количеству синхронизаций

## 8. Выводы

В ходе работы была реализована OpenMP версия алгоритма Фокса для умножения плотных матриц.

### Полученные результаты

- Реализован параллельный алгоритм умножения матриц
- Успешно пройдены все функциональные тесты
- Достигнуто ускорение более чем в 10 раз
- Реализована блочная схема вычислений

### Особенности реализации

- Использование OpenMP директив
- Блочная организация данных
- Улучшенная кэш-локальность
- Автоматическое определение размера блока

### Ограничения

- Производительность зависит от размера блока
- Возможны накладные расходы OpenMP
- Ограничение масштабирования пропускной способностью памяти

### Перспективы развития

- Автоматический подбор оптимального размера блока
- Добавление NUMA-оптимизаций
- Использование SIMD-инструкций
- Реализация гибридной MPI + OpenMP версии

## 9. Источники

1. Лекции по параллельному программированию Сысоева А. В.
2. Документация OpenMP: https://www.openmp.org/
3. Материалы курса: https://github.com/learning-process/ppc-2026-threads
4. Fox G. C. Matrix Algorithms on Parallel Hardware.
5. OpenMP Application Programming Interface Version 5.2

## 10. Приложение

```cpp
#include "sinev_a_mult_matrix_fox_algorithm/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmOMP::SinevAMultMatrixFoxAlgorithmOMP(
    const InType &in) {

  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmOMP::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 &&
         matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmOMP::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  GetOutput() =
      std::vector<double>(matrix_size * matrix_size, 0.0);

  return true;
}

void SinevAMultMatrixFoxAlgorithmOMP::SimpleMultiply(
    size_t n,
    const std::vector<double> &a,
    const std::vector<double> &b,
    std::vector<double> &c) {

#pragma omp parallel for default(none) shared(n, a, b, c) collapse(2)

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {

      double sum = 0.0;

      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] *
               b[(k * n) + j];
      }

      c[(i * n) + j] = sum;
    }
  }
}

void SinevAMultMatrixFoxAlgorithmOMP::DecomposeToBlocks(
    const std::vector<double> &src,
    std::vector<double> &dst,
    size_t n,
    size_t bs,
    int q) {

#pragma omp parallel for default(none) \
shared(src, dst, n, bs, q) collapse(2)

  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {

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
}

void SinevAMultMatrixFoxAlgorithmOMP::AssembleFromBlocks(
    const std::vector<double> &src,
    std::vector<double> &dst,
    size_t n,
    size_t bs,
    int q) {

#pragma omp parallel for default(none) \
shared(src, dst, n, bs, q) collapse(2)

  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {

      const size_t block_off =
          (static_cast<size_t>((bi * q) + bj)) *
          (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {

          const size_t src_idx =
              block_off + (i * bs) + j;

          const size_t dst_idx =
              ((static_cast<size_t>(bi) * bs + i) * n) +
              (static_cast<size_t>(bj) * bs + j);

          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmOMP::FoxStep(
    const std::vector<double> &blocks_a,
    const std::vector<double> &blocks_b,
    std::vector<double> &blocks_c,
    size_t bs,
    int q,
    int step) {

  const size_t block_size_bytes = bs * bs;

#pragma omp parallel for default(none) \
shared(blocks_a, blocks_b, blocks_c,
       bs, q, step, block_size_bytes) collapse(2)

  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) {

      const int k = (i + step) % q;

      const size_t a_off =
          (static_cast<size_t>((i * q) + k)) *
          block_size_bytes;

      const size_t b_off =
          (static_cast<size_t>((k * q) + j)) *
          block_size_bytes;

      const size_t c_off =
          (static_cast<size_t>((i * q) + j)) *
          block_size_bytes;

      for (size_t ii = 0; ii < bs; ++ii) {
        for (size_t kk = 0; kk < bs; ++kk) {

          const double val =
              blocks_a[a_off + (ii * bs) + kk];

          for (size_t jj = 0; jj < bs; ++jj) {

            blocks_c[c_off + (ii * bs) + jj] +=
                val *
                blocks_b[b_off + (kk * bs) + jj];
          }
        }
      }
    }
  }
}

bool SinevAMultMatrixFoxAlgorithmOMP::RunImpl() {
  const auto &input = GetInput();

  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  auto &c = GetOutput();

  // Для маленьких матриц используем обычное умножение
  if (n <= 8) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  size_t bs = 1;

  auto sqrt_n =
      static_cast<size_t>(
          std::sqrt(static_cast<double>(n)));

  for (size_t div = sqrt_n; div >= 1; --div) {
    if (n % div == 0) {
      bs = div;
      break;
    }
  }

  const int actual_q = static_cast<int>(n / bs);

  auto total_blocks =
      static_cast<size_t>(actual_q) *
      static_cast<size_t>(actual_q);

  auto block_elements = bs * bs;

  std::vector<double> blocks_a(
      total_blocks * block_elements);

  std::vector<double> blocks_b(
      total_blocks * block_elements);

  std::vector<double> blocks_c(
      total_blocks * block_elements,
      0.0);

  DecomposeToBlocks(
      a,
      blocks_a,
      n,
      bs,
      actual_q);

  DecomposeToBlocks(
      b,
      blocks_b,
      n,
      bs,
      actual_q);

  for (int step = 0; step < actual_q; ++step) {
    FoxStep(
        blocks_a,
        blocks_b,
        blocks_c,
        bs,
        actual_q,
        step);
  }

  AssembleFromBlocks(
      blocks_c,
      c,
      n,
      bs,
      actual_q);

  return true;
}

bool SinevAMultMatrixFoxAlgorithmOMP::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm
```