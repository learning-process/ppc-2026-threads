# Поразрядная сортировка для вещественных чисел (тип double) с простым слиянием

- **Студент:** Крымова Кристина Дмитриевна
- **Группа:** 3823Б1Пмоп3
- **Технология:** SEQ
- **Вариант:** 19

---

## 1. Контекст

Рассматривается задача поразрядной сортировки (LSD Radix Sort) для массива
вещественных чисел типа `double` с последующим слиянием.
Размер массива: **20,000,000 элементов** (160 MB).

Последовательная реализация используется как:

- эталон корректности,
- база для сравнения ускорений в параллельных версиях.

---

## 2. Постановка задачи

Дан массив `A` размера `n` с элементами типа `double`. Необходимо:

- отсортировать массив по возрастанию,
- сохранить результат в новом массиве.

### Условия корректности

- Проверяется, что массив не пуст,
- Результат должен быть отсортирован по неубыванию,
- Допустимая погрешность: **1e-14**,
- Сравнение с эталоном.

---

## 3. Базовый алгоритм

Используется двухфазный подход:

### Фаза 1. Преобразование double → uint64_t

Для корректной сортировки отрицательных чисел выполняется преобразование
IEEE 754 double в 64-битное целое с сохранением порядка.

### Фаза 2. LSD Radix Sort

Выполняется 8 проходов по 8 бит (всего 64 бита):

- подсчёт гистограммы,
- префиксные суммы,
- перестановка элементов.

### Фаза 3. Обратное преобразование uint64_t → double

### Фаза 4. Итеративное слияние

- массив разбивается на блоки фиксированного размера,
- каждый блок сортируется поразрядно,
- выполняется последовательное слияние соседних блоков.

### Асимптотика

- Время: **O(n log n)**
- Память: **O(n)**

---

## 4. Детали реализации

**Файлы:**

- `seq/include/ops_seq.hpp`
- `seq/src/ops_seq.cpp`

**Pipeline реализации:**

- **ValidationImpl** — проверка, что массив не пуст,
- **PreProcessingImpl** — копирование входных данных в выходной массив,
- **RunImpl** — вызов `IterativeMergeSort`,
- **PostProcessingImpl** — проверка корректности сортировки.

### Ключевые функции

```cpp
void LSDSortDouble(double* arr, int size);
void MergeSections(double* left, const double* right,
                   int left_size, int right_size);
void SortSections(double* arr, int size, int portion);
void IterativeMergeSort(double* arr, int size, int portion);
```

## 5. Проверка корректности

Корректность проверялась:

- сравнением с эталонной реализацией,
- функциональными тестами:
  - массивы малых размеров,
  - граничные случаи (пустой массив, один элемент),
  - случайные массивы,
  - отсортированные и обратно отсортированные массивы,
  - массивы с отрицательными числами.

**Точность сравнения:** `1e-14`

---

## 6. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| **CPU** | Apple M1 (8 ядер) |
| **RAM** | 16 GB |
| **OS** | macOS 15.6.1 |
| **Компилятор** | Clang 21.1.7 (ARM64) |
| **Сборка** | Release |

### Сборка и запуск

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cd build/bin
```

#### Функциональные тесты

```bash
./ppc_func_tests --gtest_filter="*krymova_k_lsd_sort_merge_double_seq_enabled*"
```

#### Производительность

```bash
./ppc_perf_tests --gtest_filter="*task_run_krymova_k_lsd_sort_merge_double_seq_enabled"
```

## 7. Результаты

### Параметры эксперимента

| Параметр | Значение |
| --- | --- |
| Размер массива | 20,000,000 элементов |
| Тип данных | `double` |
| Диапазон значений | случайные числа от -10000 до 10000 |

### Базовое время выполнения (baseline)

| workers | time (s) | speedup | efficiency |
| --- | --- | --- | --- |
| 1 | 0.539 | 1.00 | 1.00 |

### Наблюдения

- Основное время выполнения сосредоточено в фазе поразрядной сортировки
- Производительность ограничена пропускной способностью памяти
- Полученное время используется как базовое для сравнения
  с параллельными версиями

---

## 8. Выводы

- Реализован последовательный алгоритм поразрядной сортировки
  для типа `double`
- Корректность подтверждена функциональными тестами
- Получено базовое время выполнения **0.539 секунды**
  для массива 20M элементов
- Используемая схема (LSD Radix Sort + итеративное слияние)
  подготавливает алгоритм к эффективному распараллеливанию

---

## Список литературы

1. LSD Radix Sort:
   <https://en.wikipedia.org/wiki/Radix_sort>
2. IEEE 754 стандарт:
   <https://en.wikipedia.org/wiki/IEEE_754>
3. Google Test documentation:
   <https://google.github.io/googletest/>

---

## Приложение

### Фрагмент кода — поразрядная сортировка

```cpp
void KrymovaKLsdSortMergeDoubleSEQ::LSDSortDouble(double* arr, int size) {
    if (size <= 1) return;

    const int k_bits_per_pass = 8;
    const int k_radix = 1 << k_bits_per_pass;
    const int k_passes = 8;

    std::vector<uint64_t> ull_arr(size);
    std::vector<uint64_t> ull_tmp(size);

    for (int i = 0; i < size; ++i) {
        ull_arr[i] = DoubleToULL(arr[i]);
    }

    std::vector<unsigned int> count(k_radix, 0U);

    for (int pass = 0; pass < k_passes; ++pass) {
        int shift = pass * k_bits_per_pass;
        std::ranges::fill(count, 0U);

        for (int i = 0; i < size; ++i) {
            unsigned int digit = (ull_arr[i] >> shift) & (k_radix - 1);
            ++count[digit];
        }

        for (int i = 1; i < k_radix; ++i) {
            count[i] += count[i - 1];
        }

        for (int i = size - 1; i >= 0; --i) {
            unsigned int digit = (ull_arr[i] >> shift) & (k_radix - 1);
            ull_tmp[--count[digit]] = ull_arr[i];
        }

        ull_arr.swap(ull_tmp);
    }

    for (int i = 0; i < size; ++i) {
        arr[i] = ULLToDouble(ull_arr[i]);
    }
}
```
