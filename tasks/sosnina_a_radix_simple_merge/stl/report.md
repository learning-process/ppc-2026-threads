# Поразрядная сортировка целых чисел с простым слиянием — STL 

- **Student:** Соснина Александра Антоновна, группа 3823Б1ПР1  
- **Technology:** STL / `std::thread`  
- **Вариант:** № 17

---

## 1. Контекст

Параллелизм задаётся **явно**: шаблон **`ParallelForRange`** режет диапазон индексов `[begin, end)` на блоки, для каждого блока создаётся `std::thread`, в конце выполняется **`join`** для всех потоков. Так реализованы и фаза radix по частям, и фаза merge по парам на уровне дерева — **без** `std::execution::par` и без OpenMP/TBB.

## 2. Постановка задачи

**Цель.** Дан непустой одномерный массив целых чисел. Требуется построить **перестановку** тех же элементов, упорядоченную по **неубыванию**.

**Типы.** В `common/include/common.hpp`: **`InType` = `OutType` = `std::vector<int>`**, **`BaseTask` = `ppc::task::Task<InType, OutType>`**.

**Вход.** Вектор `data` типа `InType`, **`|data| ≥ 1`**. Пустой вход отклоняется в **`ValidationImpl()`**.

**Выход.** Вектор `OutType` той же длины, **отсортированный по неубыванию**, с совпадающим с входом **мультимножеством** значений.

**Проверка корректности.** Общие функциональные тесты в `tests/functional/main.cpp` (**48** кейсов), сравнение с эталоном **`std::sort`**; в **`RunImpl()`** дополнительно **`std::ranges::is_sorted`**.

**Особые случаи.** При **`|data| ≤ 1`** сортировка тривиальна.

Постановка **идентична** корневому `report.md` и ветке **SEQ**; в ветке **STL** меняется только способ распараллеливания (**`std::thread`**, `ParallelForRange`), а не математическая формулировка задачи.

## 3. Базовый алгоритм

Как в других ветках: radix LSD по частям, затем дерево слияний. Отличие — только способ распараллеливания циклов по индексам частей и пар.

## 4. Схема распараллеливания

1. **`num_parts`** вычисляется так же, как в TBB (пороги `kMinElementsPerPart`, `kLargeArrayThreshold`, `per_thread_floor`, `max_parts_by_grain`).
2. **Radix:** `ParallelForRange(0, num_parts, num_threads, …)` — каждый индекс `i` обрабатывает `parts[i]` с локальным буфером; разные `i` не разделяют память записи результата radix.
3. **Merge:** `ParallelForRange(0, pair_count, num_threads, …)` — для каждой пары индексов слияние в `next[idx]`; после merge исходные векторы пары освобождаются.
4. **Потоки:** `join` вызывается **после** создания всех потоков данного `ParallelForRange` (сначала цикл `emplace_back`, затем цикл `join`) — потоки реально работают параллельно внутри одного вызова `ParallelForRange`.

**Синхронизация:** между частями нет общих записей без разделения; на merge каждый `idx` пишет в свой `next[idx]`.

## 5. Детали реализации

- **Файлы:** `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`, класс `SosninaATestTaskSTL`.
- **`ParallelForRange`:** число потоков ограничивается `min(num_threads, end - begin)`.

Фрагмент запуска потоков:

```77:92:tasks/sosnina_a_radix_simple_merge/stl/src/ops_stl.cpp
  const size_t chunk = (n + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);
  std::vector<std::thread> threads;
  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const size_t lo = begin + (static_cast<size_t>(thread_idx) * chunk);
    // ...
    threads.emplace_back([lo, hi, &func]() {
      for (size_t i = lo; i < hi; ++i) {
        func(i);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
```

## 6. Проверка корректности

Те же **48** кейсов; сравнение с SEQ; `std::ranges::is_sorted` в конце `RunImpl()`.

## 7. Экспериментальная среда

- **`PPC_NUM_THREADS`**, **`OMP_NUM_THREADS`**; для раннера **`PPC_NUM_PROC=1`**.
- **Среда замеров:** Apple M2, macOS, 16 ГБ, Apple clang, Release.

Накладные расходы: создание потоков на **каждом** вызове `ParallelForRange` (в т. ч. на уровнях merge).

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
| 2 | 1,86 | 93% |
| 4 | 2,55 | 63,8% |
| 8 | 2,92 | 36,5% |

### 8.2. Режим **`pipeline`**

| N | S vs SEQ | Eff |
| - | -------- | --- |
| 2 | 1,71 | 85,5% |
| 4 | 2,34 | 58,5% |
| 8 | 2,68 | 33,5% |

Сводные таблицы и сравнение режимов — **`report.md`**, §5.

## 9. Выводы

**Потоки и синхронизация.** `ParallelForRange` создаёт пул **`std::thread`**, выполняет работу по блокам индексов и вызывает **`join`** только **после** создания всех потоков — это исключает сериализацию из-за **`join` внутри цикла** `emplace_back`. **`mutex`**, **`atomic`** и **`condition_variable`** в данной реализации **не используются**: разные потоки пишут в **разные** элементы **`next`** / в **разные** части **`parts`**, гонок по общим аккумуляторам нет.

**Накладные расходы.** На каждом уровне дерева merge вызывается **новый** `ParallelForRange`, то есть потоки **пересоздаются**; на **20 млн** элементах это всё ещё окупается (§8.1), но **`S`** чуть ниже, чем у TBB. В режиме **`pipeline`** (§8.2) дополнительно измеряется весь конвейер задачи, поэтому **`Eff`** ниже, чем в **`task_run`**.

**Итог.** STL-ветка демонстрирует **явный** контроль над параллелизмом и хорошую **переносимость**; по **`S`** она **между** OMP и TBB, что ожидаемо. Локальная схема совпадает по идее с **`SortLocalStlParallel`** в **ALL**. Итоговые рекомендации по выбору технологии — в корневом **`report.md`**.

## 10. Источники

1. Документация курса.  
2. [cppreference — std::thread](https://en.cppreference.com/w/cpp/thread/thread).  
3. [cppreference — std::ranges::merge](https://en.cppreference.com/w/cpp/algorithm/ranges/merge).
