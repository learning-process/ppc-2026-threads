# Отчет TBB: сортировка Хоара с простым слиянием

## Контекст и базовый алгоритм

TBB-версия сортирует массив `int` блоками по 64 элемента, затем сливает соседние
отсортированные диапазоны (`tbb/src/ops_tbb.cpp:18`,
`tbb/src/ops_tbb.cpp:127`). Локальная сортировка использует
ту же схему Хоара, что и `seq`: pivot из середины, два индекса и обмен до
пересечения (`tbb/src/ops_tbb.cpp:28`).

## TBB-примитивы

Используется `oneapi::tbb::parallel_for` с `oneapi::tbb::blocked_range<size_t>`
для диапазона индексов блоков и индексов слияния
(`tbb/src/ops_tbb.cpp:117`,
`tbb/src/ops_tbb.cpp:131`). Grainsize явно не передан,
поэтому применяется значение конструктора `blocked_range` по умолчанию;
partitioner также явно не указан, значит используется стандартное разбиение
`parallel_for`. Конкуренция ограничивается раннером через
`tbb::global_control(max_allowed_parallelism, ppc::util::GetNumThreads())`
(`modules/runners/src/runners.cpp:150`);
`GetNumThreads` читает `PPC_NUM_THREADS`
(`modules/util/src/util.cpp:23`).

Фрагмент, `tbb/src/ops_tbb.cpp:117`: `blocked_range` задает
независимые номера блоков.

```cpp
oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<size_t>(0, block_count),
    [this, size](const auto &range) {
  for (size_t block_index = range.begin(); block_index != range.end();
       ++block_index) {
    size_t block_start = block_index * kBlockSize;
    size_t block_end = std::min(block_start + kBlockSize, size);
    if ((block_end - block_start) > 1) {
      HoareQuickSort(data_, static_cast<int>(block_start),
                     static_cast<int>(block_end - 1));
    }
  }
});

for (size_t merge_width = kBlockSize; merge_width < size; merge_width *= 2) {
  std::vector<int> merged_data(size);
  const size_t merge_count = (size + (2 * merge_width) - 1) / (2 * merge_width);

  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, merge_count),
      [this, size, merge_width, &merged_data](const auto &range) {
```

## Детали pipeline и гонки

`ValidationImpl`, `PreProcessingImpl`, `RunImpl`, `PostProcessingImpl`
расположены в `tbb/src/ops_tbb.cpp:99`. Сортировка блоков
пишет в непересекающиеся отрезки `data_`; слияние пишет в непересекающиеся
отрезки `merged_data`, читая `data_`
(`tbb/src/ops_tbb.cpp:138`). `data_.swap` выполняется после
завершения `parallel_for` (`tbb/src/ops_tbb.cpp:148`).

## Корректность и среда

Функциональный тест сравнивает результат с `std::ranges::sort`
(`tests/functional/main.cpp:35`). Запуск
текущего `build_olesnitskiy/bin/ppc_func_tests` прошел для `seq/omp/stl/tbb`: 60
passed; ALL отдельно прошел под `mpirun -np 2`: 15 passed. Performance-вход:
`N=100000` (`tests/performance/main.cpp:20`).

## Результаты

Baseline: `seq` `TaskRun = 0.0058254364 s`; для pipeline baseline `0.0068995056
s`. Framework выполняет 5 повторов по умолчанию
(`modules/performance/include/performance.hpp:21`).

- threads: 1; time: 0.0024417256 s; speedup: 2.386; efficiency: 2.386; notes:
  `TaskRun`, `PPC_NUM_THREADS=1`.
- threads: 2; time: 0.0014976288 s; speedup: 3.889; efficiency: 1.945; notes:
  `TaskRun`, `PPC_NUM_THREADS=2`.
- threads: 4; time: 0.0011774682 s; speedup: 4.947; efficiency: 1.237; notes:
  `TaskRun`, `PPC_NUM_THREADS=4`.
- threads: 1; time: 0.0067520042 s; speedup: 1.022; efficiency: 1.022; notes:
  `pipeline`, `PPC_NUM_THREADS=1`.
- threads: 4; time: 0.0026312252 s; speedup: 2.622; efficiency: 0.656; notes:
  `pipeline`, `PPC_NUM_THREADS=4`.

## Выводы

В измерениях TBB показал лучший результат среди потоковых backend-ов:
`0.0011774682 s` при 4 потоках, speedup `4.947`. Числа относятся к `N=100000` и
фиксируются командой `ppc_perf_tests`; вывод не распространяется на другие
размеры без дополнительных замеров.
