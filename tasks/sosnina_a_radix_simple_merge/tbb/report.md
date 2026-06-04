# Поразрядная сортировка целых чисел с простым слиянием — oneTBB

- **Student:** Соснина Александра Антоновна, группа 3823Б1ПР1  
- **Technology:** oneTBB  
- **Вариант:** № 17

---

## 1. Контекст

Реализация использует **`tbb::parallel_for`** с **`tbb::simple_partitioner`** для параллельной фазы radix по частям и
для параллельного слияния пар на каждом уровне дерева. Число частей ограничивается порогами зернистости и
`GetNumThreads()`, что уменьшает избыточную глубину merge на больших входах (в т. ч. при **20 млн** элементов в
perf-тесте).

## 2. Постановка задачи

**Цель.** Дан непустой одномерный массив целых чисел. Требуется построить **перестановку** тех же элементов,
упорядоченную по **неубыванию**.

**Типы.** В `common/include/common.hpp`: **`InType` = `OutType` = `std::vector<int>`**, **`BaseTask` =
`ppc::task::Task<InType, OutType>`**.

**Вход.** Вектор `data` типа `InType`, **`|data| ≥ 1`**. Пустой вход отклоняется в **`ValidationImpl()`**.

**Выход.** Вектор `OutType` той же длины, **отсортированный по неубыванию**, с совпадающим с входом
**мультимножеством** значений.

**Проверка корректности.** Общие функциональные тесты в `tests/functional/main.cpp` (**48** кейсов), сравнение с
эталоном **`std::sort`**; в **`RunImpl()`** дополнительно **`std::ranges::is_sorted`**.

**Особые случаи.** При **`|data| ≤ 1`** сортировка тривиальна.

Постановка **идентична** корневому `report.md` и ветке **SEQ**; в ветке **TBB** изменяется только реализация
**`RunImpl`** (примитивы oneTBB и декомпозиция на части), а не формулировка задачи.

## 3. Базовый алгоритм

LSD radix (8 бит × 4 прохода, XOR для знака) по каждой части; затем итеративное попарное слияние с
`std::ranges::merge` в заранее выделенный буфер.

## 4. Схема распараллеливания

1. **Выбор `num_parts`:** базовый минимум размера части зависит от `n` и порогов `kSmallArrayThreshold` /
   `kLargeArrayThreshold`; дополнительно `per_thread_floor` (на больших `n` — порядка `n/T` потоков), затем
   `min_chunk`, `max_parts_by_grain`, и **`num_parts = min(num_threads, |data|, max_parts_by_grain)`**.
2. **Radix:** `tbb::parallel_for(0, num_parts, …, tbb::simple_partitioner{})` — каждая итерация сортирует одну часть с
   локальным буфером.
3. **Merge:** на каждом уровне `parallel_for` по `pair_count` с тем же **`simple_partitioner`**.

**`simple_partitioner`** задаёт предсказуемое разбиение диапазона итераций (удобно для воспроизводимости замеров
относительно динамического планирования по умолчанию).

## 5. Детали реализации

- **Файлы:** `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`, класс `SosninaATestTaskTBB`.

Пример вызова `parallel_for` на фазе radix:

```121:124:tasks/sosnina_a_radix_simple_merge/tbb/src/ops_tbb.cpp
  tbb::parallel_for(0, num_parts, [&](int i) {
    std::vector<int> buffer(parts[static_cast<size_t>(i)].size());
    RadixSortLSD(parts[static_cast<size_t>(i)], buffer);
  }, tbb::simple_partitioner{});
```

## 6. Проверка корректности

Общие функциональные тесты + `std::ranges::is_sorted` в `RunImpl()`.

## 7. Экспериментальная среда

- **`PPC_NUM_THREADS`**, **`OMP_NUM_THREADS`** — как в раннере; **`PPC_NUM_PROC=1`** для
  `scripts/run_tests.py --running-type=threads`.
- **Среда замеров:** Apple M2, macOS, 16 ГБ, Apple clang, Release (детали — корневой `report.md`).

```bash
export PPC_NUM_THREADS=4
export PPC_NUM_PROC=1
export OMP_NUM_THREADS=4
scripts/run_tests.py --running-type=threads
```

## 8. Результаты

### 8.1. Режим **`task_run`**

| N | Speedup vs SEQ | Efficiency |
| - | -------------- | ---------- |
| 2 | 1,91 | 95,5% |
| 4 | 2,64 | 66% |
| 8 | 3,05 | 38,1% |

### 8.2. Режим **`pipeline`**

| N | S vs SEQ | Eff |
| - | -------- | --- |
| 2 | 1,76 | 88% |
| 4 | 2,43 | 60,8% |
| 8 | 2,80 | 35% |

Корневой **`report.md`**, §5 — сводка OMP/TBB/STL для обоих режимов.

## 9. Выводы

**Примитивы и зернистость.** В коде используется **`tbb::parallel_for`** по целочисленному диапазону
**`[0, num_parts)`** с **`simple_partitioner`**, а «зернистость» выражается через расчёт **`num_parts`** (пороги
**`kSmallArrayThreshold`**, **`kLargeArrayThreshold`**, **`per_thread_floor`**, **`min_chunk`**).

**Поведение на 20 млн элементов.** При **`task_run`** (§8.1) TBB даёт **наилучшее** среди трёх потоковых веток **`S`**
при **`N = 2, 4, 8`**, что согласуется с более «толстыми» частями radix и меньшей глубиной merge относительно OMP при
больших **`T`**.

**`pipeline`.** В §8.2 **`S`** ниже, чем в **`task_run`**, за счёт полного конвейера задачи; **`Eff`** падает
пропорционально.

**Итог.** oneTBB здесь — **сильный компромисс** между выразительностью параллельного кода и контролем числа частей; при
необходимости дальнейшей настройки можно было бы экспериментировать с другим partitioner, но **`simple_partitioner`**
выбран для предсказуемости замеров. Сводные выводы по всем технологиям — в **`report.md`**.

## 10. Источники

1. Документация курса.  
2. [oneAPI Threading Building Blocks](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html).  
3. [cppreference](https://en.cppreference.com/).
