# Поразрядная сортировка для вещественных чисел (тип double) с простым слиянием

- **Студент:** Крымова Кристина Дмитриевна
- **Группа:** 3823Б1Пмоп3
- **Технология:** TBB (Intel OneTBB)
- **Вариант:** 19

---

## 1. Контекст

Рассматривается параллельная реализация поразрядной сортировки массива
вещественных чисел типа `double` с использованием библиотеки
**Intel OneTBB** [1]. TBB предоставляет высокоуровневые алгоритмы
параллельной обработки, такие как `parallel_for`, которые автоматически
распределяют работу между потоками.

**Цель:**

- реализовать параллельную версию при помощи TBB,
- исследовать масштабируемость алгоритма.

---

## 2. Постановка задачи

Дан массив `A` размера `n` с элементами типа `double`.
Необходимо отсортировать массив по возрастанию.

### Условия корректности

- Проверяется, что массив не пуст,
- Результат должен быть отсортирован по неубыванию,
- Погрешность: **1e-14**,
- Сравнение с SEQ.

---

## 3. Базовый алгоритм

Используется тот же двухфазный алгоритм, что и в SEQ:

1. Преобразование double → uint64_t,
2. LSD Radix Sort (8 проходов),
3. Обратное преобразование uint64_t → double,
4. Итеративное слияние.

**Параллелизм вводится с помощью TBB:**

- **Сортировка блоков** — `tbb::parallel_for` для распределения блоков
  между потоками,
- **Слияние** — `tbb::parallel_for` для параллельного слияния пар блоков
  на каждом уровне.

---

## 4. Схема распараллеливания

### Область параллелизма

Параллелизуются циклы по блокам (сортировка) и по парам блоков (слияние).

### Разделение данных

- `arr`, `size`, `portion` — общие (shared),
- `temp` — локальная внутри функций слияния.

### Синхронизация

- Отсутствует — TBB автоматически управляет распределением задач,
- Нет явных `mutex` или `atomic`.

### Барьеры

- Неявная синхронизация в конце `tbb::parallel_for`.

### Распределение работы

```cpp
tbb::parallel_for(tbb::blocked_range<int>(0, num_blocks, grain_size), 
    [&](const tbb::blocked_range<int> &range) {
        for (int block_index = range.begin(); block_index != range.end(); ++block_index) {
            int start = block_index * portion;
            int current_size = std::min(portion, size - start);
            LSDSortDouble(arr + start, current_size);
        }
    }
);
```

### Планировщик (partitioner)

Планировщик явно не задаётся, используется поведение TBB по умолчанию
(`auto_partitioner`), которое адаптивно распределяет итерации между
потоками.

---

## 5. Детали реализации

**Файлы:**

- `tbb/include/ops_tbb.hpp`
- `tbb/src/ops_tbb.cpp`

**Pipeline:**

- **ValidationImpl** — проверка, что массив не пуст,
- **PreProcessingImpl** — копирование входных данных,
- **RunImpl** — вызов `IterativeMergeSort` с использованием TBB,
- **PostProcessingImpl** — проверка корректности.

### Особенности

- Используется `tbb::parallel_for` для распределения блоков,
- Задан `grain_size = 1` для сортировки и `grain_size = 16` для слияния
  (экспериментально подобранные значения),
- `tbb::this_task_arena::max_concurrency()` определяет количество потоков.

---

## 6. Проверка корректности

Корректность проверялась:

- сравнением с реализацией SEQ на различных данных:
  - граничные случаи,
  - случайные массивы,
  - отсортированные и обратно отсортированные массивы.

Сравнение выполняется с точностью **1e-14**.

---

## 7. Экспериментальная среда

- **CPU**: Apple M1 (8 ядер)
- **RAM**: 16 GB
- **OS**: macOS 15.6.1
- **Компилятор**: Clang 21.1.7
- **TBB**: 2022.3.0
- **Сборка**: Release

### Сборка и запуск

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cd build/bin

export PPC_NUM_THREADS=<threads>
./ppc_perf_tests --gtest_filter="*task_run_krymova_k_lsd_sort_merge_double_tbb_enabled"
```

---

## 8. Результаты

### Параметры эксперимента

| Параметр | Значение |
| --- | --- |
| Размер массива | 20,000,000 элементов |
| Тип данных | `double` |
| Диапазон значений | случайные числа |

### Результаты измерений

| threads | time (s) | speedup | efficiency |
| --- | --- | --- | --- |
| 1 | 0.567 | 0.95 | 0.95 |
| 2 | 0.437 | 1.23 | 0.62 |
| 3 | 0.394 | 1.37 | 0.46 |
| 6 | 0.362 | 1.49 | 0.25 |
| 12 | 0.359 | 1.50 | 0.13 |

### Наблюдения

- TBB показывает **ускорение до 1.50x** на 12 потоках
- На 1 потоке TBB работает медленнее SEQ (0.95x)
  из-за накладных расходов на библиотеку
- В отличие от OMP и STL, TBB демонстрирует **устойчивую масштабируемость**
  — ускорение растёт с 1 до 12 потоков
- Однако абсолютное ускорение (1.50x) уступает OMP (1.99x) и STL (2.03x)

---

## 9. Выводы

- Реализован параллельный алгоритм поразрядной сортировки
  с использованием Intel TBB
- Корректность подтверждена функциональными тестами
- Получено ускорение до **1.50x**
- TBB показывает хорошую масштабируемость, но уступает OpenMP и std::thread
  в абсолютном ускорении
- Основным ограничением производительности является
  пропускная способность памяти

---

## Список литературы

1. Документация Intel OneTBB:
   <https://www.intel.com/content/www/us/en/docs/onetbb/get-started-guide/2023-0.html>
2. LSD Radix Sort:
   <https://en.wikipedia.org/wiki/Radix_sort>

---

## Приложение

### Фрагмент кода — параллельная сортировка блоков с TBB

```cpp
void KrymovaKLsdSortMergeDoubleTBB::SortSectionsParallel(double* arr, int size, int portion) {
    int num_blocks = (size + portion - 1) / portion;
    const size_t grain_size = 1;
    
    tbb::parallel_for(tbb::blocked_range<int>(0, num_blocks, grain_size),
        [&](const tbb::blocked_range<int> &range) {
            for (int block_index = range.begin(); block_index != range.end(); ++block_index) {
                int start = block_index * portion;
                int current_size = std::min(portion, size - start);
                LSDSortDouble(arr + start, current_size);
            }
        }
    );
}
```

### Фрагмент кода — итеративное слияние с TBB

```cpp
void KrymovaKLsdSortMergeDoubleTBB::IterativeMergeSort(double* arr, int size, int portion) {
    if (size <= 1) return;
    
    SortSectionsParallel(arr, size, portion);
    
    for (int merge_size = portion; merge_size < size; merge_size *= 2) {
        int num_pairs = (size + 2 * merge_size - 1) / (2 * merge_size);
        const size_t grain_size = 16;
        
        tbb::parallel_for(tbb::blocked_range<int>(0, num_pairs, grain_size),
            [&](const tbb::blocked_range<int> &range) {
                for (int pair_index = range.begin(); pair_index != range.end(); ++pair_index) {
                    int start = pair_index * 2 * merge_size;
                    int left_size = merge_size;
                    int right_size = std::min(merge_size, size - (start + merge_size));
                    if (right_size > 0) {
                        MergeSections(arr + start, arr + start + left_size, left_size, right_size);
                    }
                }
            }
        );
    }
}
```
