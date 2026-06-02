# Сортировка Хоара с простым слиянием - TBB

- **Студент:** Юшкова Полина Александровна, 3823Б1ПР2
- **Технология:** TBB
- **Вариант:** 13

## 1. Контекст

TBB-версия переносит две независимые фазы алгоритма на `oneapi::tbb`: сортировку блоков и слияние соседних
отсортированных диапазонов. В этой реализации я проверяю, насколько эффективно планировщик TBB распределяет работу на
входе `N=100000`.

## 2. Постановка задачи

- **Входные данные:** непустой объект `std::vector<int>`.
- **Выходные данные:** отсортированный по неубыванию `std::vector<int>`.
- **Baseline:** последовательная версия с временем `T_seq(task_run) = 0.0030497600 s`,
  `T_seq(pipeline) = 0.0083235000 s`.

## 3. Базовый алгоритм

Массив `int` сортируется блоками по 64 элемента, затем соседние отсортированные диапазоны объединяются простым слиянием.
Локальная сортировка использует схему Хоара: pivot из середины, два индекса и обмен до пересечения.

## 4. Схема распараллеливания

Используется `oneapi::tbb::parallel_for` с `oneapi::tbb::blocked_range<size_t>`:

- для сортировки блоков по индексам `block_index`;
- для каждого уровня слияния - по индексам `merge_index` независимых пар диапазонов.

Grainsize и partitioner явно не задаются, используются значения по умолчанию TBB.

Фрагмент TBB-части:

```cpp
const size_t block_count = (size + kBlockSize - 1U) / kBlockSize;
oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0U, block_count),
                          [&data, size](const oneapi::tbb::blocked_range<size_t> &range) {
  for (size_t block_index = range.begin(); block_index != range.end(); ++block_index) {
    SortBlockIfNeeded(data, size, block_index);
  }
});

for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
  MergePass(data, size, merge_width);
}
```

## 5. Детали реализации

`ValidationImpl`, `PreProcessingImpl`, `RunImpl` и `PostProcessingImpl` расположены в `tbb/src/ops_tbb.cpp`.

- Сортировка блоков пишет в непересекающиеся отрезки `data_`.
- На каждом проходе слияния создается буфер `merged_data(size)` и параллельно заполняются непересекающиеся диапазоны:
  чтение из `data_`, запись в `merged_data`.
- После завершения `parallel_for` выполняется `data_.swap(merged_data)`.

Важно: `PreProcessingImpl` копирует вход в `data_`, а `GetOutput()` заполняется только если `data_` оказался
отсортированным (`std::ranges::is_sorted(data_)`).

Файлы реализации: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.

## 6. Проверка корректности

TBB-backend подключен в общий список задач функционального теста, а эталон строится через `std::ranges::sort`.
Внутри реализации дополнительно используется проверка `std::ranges::is_sorted`.

## 7. Экспериментальная среда

- **ОС:** Windows
- **Compiler:** C++
- **Размер входных данных:** `N=100000`
- **Диапазон значений:** `[-1000000, 1000000]`
- **Baseline TaskRun:** `0.0030497600 s`
- **Baseline pipeline:** `0.0083235000 s`
- **Число повторов:** 5 по умолчанию

## 8. Результаты

Замеры выполнены на `N=100000` (5 повторов) с ограничением конкуренции через
`oneapi::tbb::global_control(max_allowed_parallelism, PPC_NUM_THREADS)`.

- threads: 1; time: 0.0076817000 s; speedup: 0.397; efficiency: 0.397; notes: `TaskRun`, `PPC_NUM_THREADS=1`.
- threads: 2; time: 0.0048485400 s; speedup: 0.629; efficiency: 0.315; notes: `TaskRun`, `PPC_NUM_THREADS=2`.
- threads: 4; time: 0.0028551400 s; speedup: 1.068; efficiency: 0.267; notes: `TaskRun`, `PPC_NUM_THREADS=4`.
- threads: 1; time: 0.0113887600 s; speedup: 0.731; efficiency: 0.731; notes: `pipeline`, `PPC_NUM_THREADS=1`.
- threads: 4; time: 0.0037098400 s; speedup: 2.244; efficiency: 0.561; notes: `pipeline`, `PPC_NUM_THREADS=4`.

## 9. Выводы

TBB-версия реализует независимую сортировку блоков и независимое слияние на каждом уровне, поэтому масштабируется лучше
простых вариантов, если объем работы достаточно большой. Для `N=100000` итог зависит от накладных расходов `parallel_for`
и выделения временных буферов на проходах слияния; для других размеров нужны отдельные замеры.
