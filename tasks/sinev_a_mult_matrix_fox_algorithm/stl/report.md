# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Фокса.

- Студент: Синев Артём Александрович, группа 3823Б1ПР2
- Технология: STL
- Вариант: 2

## 1. Введение

Умножение плотных матриц является одной из фундаментальных задач параллельных вычислений. Операция активно применяется в задачах компьютерной графики, машинного обучения, научного моделирования и численных методов.

В данной работе реализован параллельный алгоритм умножения матриц с использованием стандартной библиотеки потоков C++ (std::thread). В основе реализации лежит блочная схема алгоритма Фокса, позволяющая эффективно распределять вычисления между потоками и улучшать локальность данных.

Основной целью работы является реализация масштабируемого алгоритма умножения матриц с использованием нативных потоков C++ и анализ эффективности по сравнению с Intel oneTBB и OpenMP.

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

В реализации используется блочная схема алгоритма Фокса с применением пула потоков std::thread.

### Основная идея

Матрицы разбиваются на квадратные блоки размера:

```math
bs \times bs
```

Каждый поток обрабатывает отдельные блоки матриц. Для распределения работы используется std::atomic счетчик для организации динамического распределения блоков.

### Выбор размера блока

```cpp
size_t ChooseBlockSize(size_t n) {
  if (n % 128 == 0) return 128;
  if (n % 64 == 0) return 64;
  if (n % 32 == 0) return 32;
  if (n % 16 == 0) return 16;
  return 1;
}
```

### Распределение блоков

```cpp
unsigned int num_threads = std::thread::hardware_concurrency();
std::atomic<size_t> next_block(0);

for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
  threads.emplace_back([&]() {
    size_t block_idx = 0;
    while ((block_idx = next_block.fetch_add(1)) < total_blocks) {
      // обработка блока block_idx
    }
  });
}
```

### Используемые механизмы 

- `std::thread` — создание и управление потоками
- `std::atomic` — атомарные операции для распределения задач
- `std::thread::hardware_concurrency()` — определение количества ядер
- Ручная балансировка нагрузки через динамическое распределение

## 5. Детали реализации

### Структура проекта

- `common/include/common.hpp`
- `stl/include/ops_stl.hpp`
- `stl/src/ops_stl.cpp`
- `tests/functional/main.cpp`
- `tests/performance/main.cpp`

### Основной класс

```cpp
class SinevAMultMatrixFoxAlgorithmSTL : public BaseTask
```

### Основные методы

- `ValidationImpl()`
- `PreProcessingImpl()`
- `RunImpl()`
- `PostProcessingImpl()`

### Проверка входных данных

```cpp
bool SinevAMultMatrixFoxAlgorithmSTL::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  return matrix_size > 0 &&
         matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}
```

### Разбиение матриц на блоки

```cpp
void SinevAMultMatrixFoxAlgorithmSTL::DecomposeToBlocks(
    const std::vector<double> &src,
    std::vector<double> &dst,
    size_t n, size_t bs, int q) {
  
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) num_threads = 2;

  std::vector<std::thread> threads;
  std::atomic<size_t> next_block(0);
  size_t total_blocks = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&]() {
      size_t block_idx = 0;
      while ((block_idx = next_block.fetch_add(1)) < total_blocks) {
        int bi = static_cast<int>(block_idx / q);
        int bj = static_cast<int>(block_idx % q);
        const size_t block_off = block_idx * (bs * bs);
        
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            size_t src_idx = ((bi * bs + i) * n) + (bj * bs + j);
            size_t dst_idx = block_off + (i * bs) + j;
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    });
  }

  for (auto &thread : threads) thread.join();
}
```

### Умножение блоков

```cpp
void SinevAMultMatrixFoxAlgorithmSTL::MultiplyBlocks(
    const std::vector<double> &blocks_a,
    const std::vector<double> &blocks_b,
    std::vector<double> &blocks_c,
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
```

### Выполнение шага алгоритма Фокса

```cpp
void SinevAMultMatrixFoxAlgorithmSTL::FoxStep(
    const std::vector<double> &blocks_a,
    const std::vector<double> &blocks_b,
    std::vector<double> &blocks_c,
    size_t bs, int q, int step) {
  
  const size_t block_size = bs * bs;
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) num_threads = 2;

  std::vector<std::thread> threads;
  std::atomic<size_t> next_cell(0);
  size_t total_cells = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&]() {
      size_t cell_idx = 0;
      while ((cell_idx = next_cell.fetch_add(1)) < total_cells) {
        int i = static_cast<int>(cell_idx / q);
        int j = static_cast<int>(cell_idx % q);
        const int k = (i + step) % q;

        size_t a_off = (static_cast<size_t>((i * q) + k)) * block_size;
        size_t b_off = (static_cast<size_t>((k * q) + j)) * block_size;
        size_t c_off = (static_cast<size_t>((i * q) + j)) * block_size;

        MultiplyBlocks(blocks_a, blocks_b, blocks_c, bs, a_off, b_off, c_off);
      }
    });
  }

  for (auto &thread : threads) thread.join();
}
```

### Основной вычислительный метод 

```cpp
bool SinevAMultMatrixFoxAlgorithmSTL::RunImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);
  auto &c = GetOutput();

  // Для маленьких матриц используем простое умножение
  if (n <= 64) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  size_t bs = 64;
  while (n % bs != 0 && bs > 16) {
    bs /= 2;
  }

  if (n % bs != 0) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

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
```

### Особенности реализации

- Используется `std::thread::hardware_concurrency()` для определения числа потоков
- Применяется атомарный счетчик для динамического распределения работы
- Реализована адаптивная логика выбора размера блока
- Для малых матриц (≤64) используется последовательное умножение
- Все потоки синхронизируются через `join()`

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
./build/bin/ppc_func_tests --gtest_filter="*Sinev*stl*" -v
```

Запуск performance тестов:

```bash
./build/bin/ppc_perf_tests --gtest_filter="*Sinev*stl*" -v
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

| Режим              | Потоки | Время, с   | Ускорение | Эффективность |
|--------------------|---------|-------------|------------|----------------|
| seq                | 1       | 0.0261574700 | 1.00       | 100%           |
| omp (pipeline)     | 16      | 0.0024294812 | 10.77      | 67.3%          |
| omp (task_run)     | 16      | 0.0042057726 | 6.22       | 38.9%          |
| tbb (pipeline)     | 16      | 0.0166006096 | 1.58       | 9.9%           |
| tbb (task_run)     | 16      | 0.0166640150 | 1.57       | 9.8%           |
| stl (pipeline)     | 16      | 0.0124331014 | 2.10       | 13.1%           |
| stl (task_run)     | 16      | 0.0104718960 | 2.50       | 15.6%          |

Ускорение:

```math
Speedup = \frac{T_{seq}}{T_{omp}}
```

Эффективность:

```math
Efficiency = \frac{Speedup}{Threads} \times 100\%
```
#### Анализ результатов

STL (std::thread):

- Реализация STL показала ускорение до 2.5 раз при использовании 16 потоков
- STL оказалась на 37-58% быстрее чем TBB
- Эффективность использования потоков составила 13-16%

### Основные причины различий

1. OpenMP лидирует благодаря:
- Низким накладным расходам на создание задач
- Эффективной работе с кэш-памятью
- Оптимизированному планировщику

2. STL занимает промежуточное положение:
- Нативные потоки имеют меньшие накладные расходы, чем TBB
- Атомарные операции обеспечивают балансировку
- Ручное управление дает больше контроля
- Нет автоматической оптимизации как в OpenMP

3. TBB показывает наихудшие результаты:
- Высокие накладные расходы task scheduler
- Дополнительная абстракция задач
- Неоптимальный для данного размера матриц

### Узкие места алгоритма
- Ручное управление потоками требует больше кода
- Нет встроенной балансировки нагрузки (реализована вручную)
- Нет автоматической настройки под архитектуру

### Масштабируемость

STL обеспечивает:
- Прямой контроль над потоками
- Минимальные накладные расходы на синхронизацию (только атомарные счетчики)
- Простота отладки и профилирования
- Отсутствие зависимостей от внешних библиотек

## 8. Выводы

В ходе работы была реализована STL версия алгоритма Фокса для умножения плотных матриц с использованием нативных потоков C++.

### Полученные результаты

- Реализован параллельный алгоритм на основе std::thread
- Успешно пройдены все функциональные тесты (13 тестов)
- Получено ускорение до 2.5 раз относительно последовательной версии
- Реализована блочная схема вычислений с динамическим распределением блоков
- STL оказалась эффективнее TBB (на 59%), но уступает OpenMP (на 148%)

### Особенности реализации

- Использование `std::atomic` для динамического распределения задач
- Оптимальный размер блока подбирается адаптивно
- Полное ручное управление жизненным циклом потоков
- Автоматическое определение числа ядер системы

### Ограничения

- Производительность зависит от количества ядер
- Возможна неравномерная загрузка при неудачном выборе блока
- Отсутствие встроенных механизмов оптимизации кэша
- Нет автоматического управления affinity

### Перспективы развития

- Использование `std::jthread` (C++20) для автоматического управления
- Реализация affinity через `std::thread::native_handle()`
- Настройка размера блока под кэш-иерархию
- Использование SIMD через интринсики компилятора
- Создание пула потоков с возможность переиспользования

## 9. Источники

1. Лекции по параллельному программированию Сысоева А. В.
2. Документация C++ Thread Support Library: https://en.cppreference.com/w/cpp/thread
3. Материалы курса: https://github.com/learning-process/ppc-2026-threads
4. Fox G. C. Matrix Algorithms on Parallel Hardware
5. Williams A. C++ Concurrency in Action, 2nd Edition

## 10. Приложение

```cpp
#include "sinev_a_mult_matrix_fox_algorithm/stl/include/ops_stl.hpp"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmSTL::SinevAMultMatrixFoxAlgorithmSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmSTL::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmSTL::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void SinevAMultMatrixFoxAlgorithmSTL::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double tmp = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += tmp * b[(k * n) + j];
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::atomic<size_t> next_block(0);
  size_t total_blocks = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&]() {
      size_t block_idx = 0;
      while ((block_idx = next_block.fetch_add(1)) < total_blocks) {
        int bi = static_cast<int>(block_idx / q);
        int bj = static_cast<int>(block_idx % q);

        const size_t block_off = block_idx * (bs * bs);
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
            const size_t dst_idx = block_off + (i * bs) + j;
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::atomic<size_t> next_block(0);
  size_t total_blocks = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&]() {
      size_t block_idx = 0;
      while ((block_idx = next_block.fetch_add(1)) < total_blocks) {
        int bi = static_cast<int>(block_idx / q);
        int bj = static_cast<int>(block_idx % q);

        const size_t block_off = block_idx * (bs * bs);
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            const size_t src_idx = block_off + (i * bs) + j;
            const size_t dst_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::MultiplyBlocks(const std::vector<double> &blocks_a,
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

void SinevAMultMatrixFoxAlgorithmSTL::FoxStep(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                                              std::vector<double> &blocks_c, size_t bs, int q, int step) {
  const size_t block_size = bs * bs;
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::atomic<size_t> next_cell(0);
  size_t total_cells = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&]() {
      size_t cell_idx = 0;
      while ((cell_idx = next_cell.fetch_add(1)) < total_cells) {
        int i = static_cast<int>(cell_idx / q);
        int j = static_cast<int>(cell_idx % q);
        const int k = (i + step) % q;

        const size_t a_off = (static_cast<size_t>((i * q) + k)) * block_size;
        const size_t b_off = (static_cast<size_t>((k * q) + j)) * block_size;
        const size_t c_off = (static_cast<size_t>((i * q) + j)) * block_size;

        MultiplyBlocks(blocks_a, blocks_b, blocks_c, bs, a_off, b_off, c_off);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

bool SinevAMultMatrixFoxAlgorithmSTL::RunImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);
  auto &c = GetOutput();

  if (n <= 64) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  size_t bs = 64;
  while (n % bs != 0 && bs > 16) {
    bs /= 2;
  }

  if (n % bs != 0) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

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

size_t SinevAMultMatrixFoxAlgorithmSTL::ChooseBlockSize(size_t n) {
  if (n % 128 == 0) {
    return 128;
  }
  if (n % 64 == 0) {
    return 64;
  }
  if (n % 32 == 0) {
    return 32;
  }
  if (n % 16 == 0) {
    return 16;
  }
  return 1;
}

bool SinevAMultMatrixFoxAlgorithmSTL::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm

```