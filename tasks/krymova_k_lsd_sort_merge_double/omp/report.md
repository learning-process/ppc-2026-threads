# Поразрядная сортировка для вещественных чисел (тип double) с простым слиянием

- **Студент:** Крымова Кристина Дмитриевна
- **Группа:** 3823Б1Пмоп3
- **Технология:** OMP
- **Вариант:** 19

---

## 1. Контекст

Рассматривается параллельная реализация поразрядной сортировки
(LSD Radix Sort)
для массива вещественных чисел типа `double` с использованием **OpenMP** [1].

Массив имеет размер 20,000,000 элементов (160 MB).

**Цель:**

- уменьшить время выполнения по сравнению с SEQ,
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

---

## 4. Схема распараллеливания

### Область параллелизма

Параллелизуется:

- цикл по блокам при сортировке,
- цикл по парам блоков при слиянии.

### Разделение данных

- `arr`, `size`, `portion`, `merge_size` — **shared** (общие для всех потоков),
- `current_size`, `left_size`, `right_size`, `pair_start` — **private**
  (локальные для каждого потока).

### Синхронизация

- Отсутствует — каждый поток работает со своей областью памяти,
- Нет `critical`, `atomic`, `reduction`.

### Барьеры

- Неявный барьер в конце `#pragma omp parallel for`.

### Распределение работы

Используется статическое распределение через `#pragma omp parallel for`
с автоматическим делением диапазонов компилятором.

```cpp
#pragma omp parallel for default(none) shared(arr, size, portion)
for (int block_start = 0; block_start < size; block_start += portion) {
    int current_size = std::min(portion, size - block_start);
    LSDSortDouble(arr + block_start, current_size);
}
```

---

## 5. Детали реализации

**Файлы:**

- `omp/include/ops_omp.hpp`
- `omp/src/ops_omp.cpp`

**Pipeline:**

- **ValidationImpl** — проверка, что массив не пуст,
- **PreProcessingImpl** — копирование входных данных,
- **RunImpl** — параллельный вызов `IterativeMergeSort`,
- **PostProcessingImpl** — проверка корректности.

### Особенности

- Используется `#pragma omp parallel for` для распределения блоков,
- Каждый поток работает со своим диапазоном памяти,
- Нет гонок данных благодаря независимым областям записи.

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
- **Компилятор**: Clang 21.1.7 с поддержкой OpenMP (`_OPENMP=202011`)
- **Сборка**: Release

### Сборка и запуск

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cd build/bin

export PPC_NUM_THREADS=<threads>
./ppc_perf_tests --gtest_filter="*task_run_krymova_k_lsd_sort_merge_double_omp_enabled"
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
| 1 | 0.273 | 1.97 | 1.97 |
| 2 | 0.271 | 1.99 | 0.99 |
| 3 | 0.272 | 1.98 | 0.66 |
| 6 | 0.273 | 1.97 | 0.33 |
| 12 | 0.272 | 1.98 | 0.17 |

### Наблюдения

- Ускорение достигает **1.99x** уже на 2 потоках
- Дальнейшее увеличение потоков не даёт прироста
  из-за memory-bound характера задачи
- Эффективность падает с ростом числа потоков
  из-за overhead и конкуренции за шину памяти

---

## 9. Выводы

- Реализован параллельный алгоритм поразрядной сортировки
  с использованием OpenMP
- Корректность подтверждена функциональными тестами
- Получено ускорение до **1.99x**
- Алгоритм не масштабируется дальше 2 потоков
  из-за ограничений пропускной способности памяти
- Основным ограничением производительности является память,
  а не вычислительные ресурсы

---

## Список литературы

1. Документация OpenMP:
   <https://www.openmp.org/>
2. LSD Radix Sort:
   <https://en.wikipedia.org/wiki/Radix_sort>

---

## Приложение

### Фрагмент кода — параллельная сортировка блоков

```cpp
void KrymovaKLsdSortMergeDoubleOMP::SortSectionsParallel(double* arr, int size, int portion) {
    #pragma omp parallel for default(none) shared(arr, size, portion)
    for (int block_start = 0; block_start < size; block_start += portion) {
        int current_size = std::min(portion, size - block_start);
        LSDSortDouble(arr + block_start, current_size);
    }
}
```

### Фрагмент кода — итеративное слияние

```cpp
void KrymovaKLsdSortMergeDoubleOMP::IterativeMergeSort(double* arr, int size, int portion) {
    if (size <= 1) return;
    
    SortSectionsParallel(arr, size, portion);
    
    for (int merge_size = portion; merge_size < size; merge_size *= 2) {
        #pragma omp parallel for default(none) shared(arr, size, merge_size)
        for (int pair_start = 0; pair_start < size; pair_start += 2 * merge_size) {
            int left_size = merge_size;
            int right_size = std::min(merge_size, size - (pair_start + merge_size));
            if (right_size > 0) {
                MergeSections(arr + pair_start, arr + pair_start + left_size, left_size, right_size);
            }
        }
    }
}
```
