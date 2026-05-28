# Поразрядная сортировка для вещественных чисел (тип double) с простым слиянием

- **Студент:** Крымова Кристина Дмитриевна
- **Группа:** 3823Б1Пмоп3
- **Технология:** STL (std::thread)
- **Вариант:** 19

---

## 1. Контекст

Рассматривается параллельная реализация поразрядной сортировки массива
вещественных чисел типа `double` с использованием стандартной библиотеки
потоков **std::thread**. Данная реализация выполняется в качестве учебной
для проверки способностей ручного распределения данных и функций потокам
без использования примитивов OpenMP и OneTBB.

**Цель:**

- реализовать параллельную версию при помощи `std::thread`,
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

**Параллелизм вводится на двух уровнях:**

- **Сортировка блоков** — каждый поток сортирует свой блок независимо,
- **Слияние** — параллельное слияние пар блоков на каждом уровне.

По сравнению с SEQ, STL-версия меняет только организацию вычислений
из-за необходимости разбить диапазон блоков на области,
обрабатываемые каждым потоком.

---

## 4. Схема распараллеливания

### Область параллелизма

1. Вычисляется количество блоков
2. Блоки распределяются между потоками
3. Каждый поток сортирует свои блоки последовательно
4. После завершения всех потоков выполняется итеративное слияние

### Разделение данных

- `arr`, `size`, `portion` — общие, передаются по ссылке,
- `temp` — локальная внутри функций слияния.

### Синхронизация

- Отсутствует — каждый поток пишет в свою область памяти,
- Нет `mutex`, `atomic`, `condition_variable`.

### Барьеры

- Неявный барьер в виде `th.join()` после заполнения контейнера потоков.

### Распределение работы

Используется ручное разбиение блоков между потоками:

```cpp
int num_blocks = (size + portion - 1) / portion;
int num_threads_used = std::min(num_threads, num_blocks);
int blocks_per_thread = (num_blocks + num_threads_used - 1) / num_threads_used;

for (int thread_id = 0; thread_id < num_threads_used; ++thread_id) {
    int start_block = thread_id * blocks_per_thread;
    int end_block = std::min(start_block + blocks_per_thread, num_blocks);
    threads.emplace_back([&, start_block, end_block]() {
        for (int block = start_block; block < end_block; ++block) {
            int start_pos = block * portion;
            int current_size = std::min(portion, size - start_pos);
            LSDSortDoubleSequential(arr + start_pos, current_size);
        }
    });
}

for (auto &th : threads) {
    th.join();
}
```

---

## 5. Детали реализации

**Файлы:**

- `stl/include/ops_stl.hpp`
- `stl/src/ops_stl.cpp`

**Pipeline:**

- **ValidationImpl** — проверка, что массив не пуст,
- **PreProcessingImpl** — копирование входных данных,
- **RunImpl** — создание потоков, назначение им функций,
  параллельная сортировка и слияние,
- **PostProcessingImpl** — проверка корректности.

### Особенности

- Используется двухфазный подход, за счёт которого избегаются гонки данных,
- Каждый поток работает со своим набором блоков,
- Слияние также выполняется параллельно на каждом уровне.

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
- **Сборка**: Release

### Сборка и запуск

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cd build/bin

export PPC_NUM_THREADS=<threads>
./ppc_perf_tests --gtest_filter="*task_run_krymova_k_lsd_sort_merge_double_stl_enabled"
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
| 1 | 0.267 | 2.02 | 2.02 |
| 2 | 0.268 | 2.01 | 1.01 |
| 3 | 0.270 | 2.00 | 0.67 |
| 6 | 0.266 | 2.03 | 0.34 |
| 12 | 0.266 | 2.03 | 0.17 |

### Наблюдения

- Ускорение достигает **2.03x** — лучший результат среди всех технологий
- Суперлинейное ускорение (2.02x на 1 потоке) объясняется лучшим
  использованием кэш-памяти
- Дальнейшее увеличение потоков не даёт прироста
  из-за memory-bound характера задачи
- Эффективность падает с ростом числа потоков

---

## 9. Выводы

- Реализован параллельный алгоритм поразрядной сортировки
  с использованием `std::thread`
- Корректность подтверждена функциональными тестами
- Получено максимальное ускорение **2.03x** — лучшее среди всех реализаций
- Алгоритм хорошо масштабируется до 2 потоков, но не более
- Основным ограничением производительности является
  пропускная способность памяти

---

## Список литературы

1. Документация std::thread:
   <https://en.cppreference.com/w/cpp/thread/thread>
2. LSD Radix Sort:
   <https://en.wikipedia.org/wiki/Radix_sort>

---

## Приложение

### Фрагмент кода — параллельная сортировка блоков с ручным распределением

```cpp
void KrymovaKLsdSortMergeDoubleSTL::SortSectionsParallel(double* arr, int size, int portion, int num_threads) {
    int num_blocks = (size + portion - 1) / portion;
    int num_threads_used = std::min(num_threads, num_blocks);
    
    if (num_threads_used <= 1) {
        for (int start = 0; start < size; start += portion) {
            int current_size = std::min(portion, size - start);
            LSDSortDoubleSequential(arr + start, current_size);
        }
        return;
    }
    
    std::vector<std::thread> threads;
    int blocks_per_thread = (num_blocks + num_threads_used - 1) / num_threads_used;
    
    for (int thread_id = 0; thread_id < num_threads_used; ++thread_id) {
        int start_block = thread_id * blocks_per_thread;
        int end_block = std::min(start_block + blocks_per_thread, num_blocks);
        threads.emplace_back([&, start_block, end_block]() {
            for (int block = start_block; block < end_block; ++block) {
                int start_pos = block * portion;
                int current_size = std::min(portion, size - start_pos);
                LSDSortDoubleSequential(arr + start_pos, current_size);
            }
        });
    }
    
    for (auto &th : threads) {
        th.join();
    }
}
```

### Фрагмент кода — итеративное слияние

```cpp
void KrymovaKLsdSortMergeDoubleSTL::IterativeMergeSort(double* arr, int size, int portion, int num_threads) {
    if (size <= 1) return;
    
    SortSectionsParallel(arr, size, portion, num_threads);
    
    for (int merge_size = portion; merge_size < size; merge_size *= 2) {
        for (int start_pos = 0; start_pos < size; start_pos += 2 * merge_size) {
            int left_size = merge_size;
            int right_size = std::min(merge_size, size - (start_pos + merge_size));
            if (right_size <= 0) continue;
            
            MergeSections(arr + start_pos, arr + start_pos + left_size, left_size, right_size);
        }
    }
}
```
