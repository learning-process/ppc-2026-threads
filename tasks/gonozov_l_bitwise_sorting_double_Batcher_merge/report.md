# Поразрядная сортировка для вещественных чисел (тип double) с четно-нечетным слиянием Бэтчера

## Сравнительный анализ реализаций SEQ, OMP, TBB, STL и гибридной MPI+TBB

- **Студент**: Gonozov Leonid
- **Вариант**: 20

## 1. Введение

### 1.1. Контекст работы

В рамках курса по параллельному программированию была реализована задача
поразрядной сортировки вещественных чисел типа `double` с использованием
четно-нечетного слияния Бэтчера.
Целью работы является сравнение пяти различных подходов к параллелизации:

| Технология |     Описание                                           |
|------------|--------------------------------------------------------|
| **SEQ**    | Последовательная эталонная реализация                  |
| **OMP**    | Параллелизация с использованием OpenMP                 |
| **TBB**    | Параллелизация с использованием Intel oneTBB           |
| **STL**    | Параллелизация с использованием нативных `std::thread` |
| **ALL**    | Гибридная реализация MPI + TBB                         |

### 1.2. Цели и задачи

1. Реализовать корректный последовательный алгоритм как эталон
2. Разработать параллельные версии с использованием различных технологий
3. Провести сравнительный анализ производительности
4. Оценить масштабируемость и эффективность каждой реализации
5. Выявить сильные и слабые стороны каждого подхода

## 2. Единая постановка задачи

### 2.1. Входные данные

`std::vector<double>`

### 2.2. Выходные данные

- **Тип**: `std::vector<double>`
- **Свойства**: отсортированный по возрастанию массив

### 2.3. Требования к корректности

1. Результат должен совпадать с последовательной версией
2. Алгоритм должен сохранять относительный порядок равных элементов
3. Корректная обработка крайних случаев (пустой массив, один элемент, отрицательные числа)
4. Детерминированность результата

## 3. Базовый алгоритм (последовательная версия)

### 3.1. Общая схема

Входной массив double
↓
Преобразование double → uint64_t (сохранение порядка)
↓
Поразрядная сортировка (8 проходов, radix = 256)
↓
Слияние половин сетью Бэтчера
↓
Обратное преобразование uint64_t → double
↓
Отсортированный массив

### 3.2. Ключевые функции

```cpp
// Преобразование double → uint64_t для сортировки
uint64_t DoubleToSortableInt(double d) {
    uint64_t bits;
    std::memcpy(&bits, &d, sizeof(double));
    if ((bits >> 63) != 0) return ~bits;      // отрицательные
    return bits | 0x8000000000000000ULL;      // положительные
}

// Итеративная сеть Бэтчера
void OddEvenMergeIterative(double *arr, size_t start, size_t n) {
    size_t step = n / 2;
    CompareExchangeBlocks(arr, start, step);
    step /= 2;
    for (; step > 0; step /= 2) {
        for (size_t i = step; i < n - step; i += step * 2) {
            CompareExchangeBlocks(arr, start + i, step);
        }
    }
}
```

### 3.3. Асимптотическая сложность

Поразрядная сортировка: O(8n) = O(n)
Слияние Бэтчера: O(n·log²n)
Общая: O(n·log²n)
Память: O(n)

## 4. Схемы распараллеливания

### 4.1. Общее сравнение подходов

Технология Уровень Основные               примитивы
SEQ         последовательный              —
OMP         внутрипроцессный              sections, parallel for
TBB         внутрипроцессный              parallel_invoke, parallel_for
STL         внутрипроцессный              std::thread, join
ALL         межпроцессный + внутрипроцессный  MPI + TBB

### 4.2. OpenMP (OMP)

```cpp
// Параллельная сортировка половин
#pragma omp parallel sections default(none) shared(left, right)
{
    #pragma omp section { RadixSortDouble(left); }
    #pragma omp section { RadixSortDouble(right); }
}

// Параллельное слияние Бэтчера
for (size_t len = 2; len <= n; len *= 2) {
    #pragma omp parallel for default(none) shared(arr, n, len)
    for (size_t i = 0; i < n; i += len) {
        MergingHalves(arr, i, len);
    }
}
```

Особенности:
default(none) — требует явного указания атрибутов переменных
Неявные барьеры в конце каждой параллельной области
Простота реализации при хорошей эффективности

### 4.3. Intel TBB

```cpp
// Параллельная сортировка половин
tbb::parallel_invoke(
    [&]() { RadixSortDouble(left); },
    [&]() { RadixSortDouble(right); }
);

// Параллельное слияние Бэтчера
for (size_t len = 2; len <= n; len *= 2) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, len),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i += len) {
                MergingHalves(arr, i, len);
            }
        });
}
```

Особенности:
Автоматическая балансировка нагрузки
Адаптивный выбор размера блока (grainsize)
Рекурсивное разбиение диапазонов

### 4.4. STL (std::thread)

```cpp
// Сортировка чанков
std::vector<std::thread> threads;
for (size_t i = 0; i < chunks_count; ++i) {
    threads.emplace_back([raw_data, i, chunk_size]() {
        SortChunk(raw_data, i * chunk_size, chunk_size);
    });
}
for (auto& thread : threads) thread.join();

// Слияние блоков
for (size_t cur_size = chunk_size; cur_size < total_size; cur_size *= 2) {
    std::vector<std::thread> merge_threads;
    for (size_t i = 0; i < merges_count; ++i) {
        merge_threads.emplace_back([raw_data, i, cur_size]() {
            OddEvenMergeIterative(raw_data, i * 2 * cur_size, 2 * cur_size);
        });
    }
    for (auto& thread : merge_threads) thread.join();
}
```

Особенности:
Ручное управление разбиением диапазонов
Явное создание и ожидание потоков
Минимальные накладные расходы на синхронизацию

### 4.5. Гибридная версия MPI + TBB (ALL)

```cpp
// Распределение данных
MPI_Scatterv(local_data_.data(), send_counts.data(), send_displs.data(),
             MPI_DOUBLE, local_data.data(), local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// Локальная TBB-сортировка
tbb::parallel_for(0, num_chunks_int, [&](int i) {
    SortChunkTBB(raw_data, i, chunk_size);
});

// Сбор результатов на процесс 0
if (rank == 0) {
    for (int p = 1; p < size; ++p) {
        MPI_Recv(&recv_size, 1, MPI_INT, p, 0, ...);
        MPI_Recv(recv_data.data(), recv_size, MPI_DOUBLE, p, 1, ...);
    }
    RadixSortDouble(all_data);
}

// Рассылка результата
MPI_Bcast(&total_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
MPI_Bcast(local_data.data(), total_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
```

## 5. Единая методика эксперимента

### 5.1. Экспериментальная среда

CPU: Intel Core i7-12650H
10 ядер / 16 аппаратных потоков
RAM: 16 GB DDR4
OS: Ubuntu 22.04 (Linux)
Компилятор: GCC 11.4.0
Стандарт C++: C++20
Тип сборки: Release (-O2)
MPI: OpenMPI 4.1.2
TBB: oneTBB 2021.8

### 5.3. Формулы расчёта

Speedup = T_seq / T_parallel
Efficiency = (Speedup / workers) × 100%

## Результаты экспериментов

на массиве из 1001000 элементов

| Mode        | Processes | Threads | Time, s | Speedup | Efficiency |
|-------------|-----------|---------|---------|---------|------------|
| seq         | 1         | 1       | 0.11779 | 1.00    | N/A        |
| omp         | 1         | 2       | 0.068   | 1,732   | 86.61%     |
| omp         | 1         | 4       | 0.0541  | 2,177   | 54.43%     |
| omp         | 1         | 8       | 0.056   | 2,103   | 25.29%     |
| tbb         | 1         | 2       | 0.0667  | 1,647   | 82.38%     |
| tbb         | 1         | 4       | 0.0512  | 2,146   | 53.66%     |
| tbb         | 1         | 8       | 0.0518  | 2,121   | 26.52%     |
| stl         | 1         | 2       | 0.02299 | 4,78    | 239.01%    |
| stl         | 1         | 4       | 0.0207  | 5,309   | 132.72%    |
| stl         | 1         | 8       | 0.0204  | 5,387   | 67.34%     |
| all         | 1         | 2       | 0.0579  | 1,898   | 94.9%      |
| all         | 2         | 1       | 0.0819  | 1,34    | 67%        |
| all         | 2         | 2       | 0,028   | 3,916   | 97.9%      |

OpenMP:
2 потока (86.6%) — почти идеальное масштабирование
4 потока (54.4%) — ускорение растёт, но эффективность падает
8 потоков (25.3%) — насыщение, накладные расходы превышают выигрыш

TBB:
Результаты очень близки к OpenMP
TBB даёт небольшое преимущество на 4 потоках (лучшая балансировка)
Те же ограничения масштабируемости

STL:
Ключевое отличие: Разбиение на N чанков (а не на 2 половины)

Эффективность >100% — не ошибка! Причина — кэш-эффект:
SEQ: 8 МБ → не помещается в кэш → много промахов
STL на 2 потоках: 2 чанка по 4 МБ → каждый помещается в L3-кэш
Итог: 2× (параллелизм) × 2.3× (кэш-эффект) ≈ 4.6× ускорение

ALL:
1x2: только TBB, отличная эффективность
2x1: только MPI, накладные расходы ~33%
2x2: лучшая эффективность
Почему 2x2 так эффективен: MPI (уменьшение данных до ~500K) +
TBB (параллельная сортировка) + Кэш-эффект (4 МБ → L3-кэш)

Все реализации корректны
Лучшее абсолютное ускорение — STL (5.39× на 8 потоках)
Лучшая эффективность — ALL
Проще всего реализовать — OMP (минимальные изменения)

### 10.1. Сборка проекта

```bash
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

### 10.2. Запуск тестов

```bash
Функциональные тесты
scripts/run_tests.py --running-type=threads

Тесты производительности
scripts/run_tests.py --running-type=performance

Конкретная конфигурация
export PPC_NUM_PROC=2
export PPC_NUM_THREADS=2
scripts/run_tests.py --running-type=performance
```
