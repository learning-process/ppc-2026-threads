# Сортировка Хоара с простым слиянием — TBB

- **Студент:** Олесницкий Владимир Тарасович, 3823Б1ПР2
- **Технология:** TBB
- **Вариант:** 13

## 1. Контекст

TBB-версия переносит две независимые фазы алгоритма на `oneapi::tbb`:
сортировку блоков и слияние соседних отсортированных диапазонов. Эта реализация
проверяет, насколько эффективно планировщик TBB распределяет работу для входа
`N=100000`.

## 2. Постановка задачи

- **Входные данные:** непустой объект `std::vector<int>`.
- **Выходные данные:** отсортированный по неубыванию `std::vector<int>`.
- **Baseline:** последовательная версия с временем `T_seq = 0.0058254364 s`.

## 3. Базовый алгоритм

Массив `int` сортируется блоками по 64 элемента, затем соседние отсортированные
диапазоны объединяются простым слиянием. Локальная сортировка использует ту же
схему Хоара, что и SEQ: pivot из середины, два индекса и обмен до пересечения.

## 4. Схема распараллеливания

Используется `oneapi::tbb::parallel_for` с
`oneapi::tbb::blocked_range<size_t>` для диапазона номеров блоков и диапазона
номеров слияний. Grainsize и partitioner явно не задаются, поэтому применяются
значения по умолчанию.

Конкуренция ограничивается раннером через
`tbb::global_control(max_allowed_parallelism, ppc::util::GetNumThreads())`.
`GetNumThreads` читает `PPC_NUM_THREADS`.

Фрагмент TBB-части:

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

## 5. Детали реализации

`ValidationImpl`, `PreProcessingImpl`, `RunImpl` и `PostProcessingImpl`
расположены в `tbb/src/ops_tbb.cpp`. Сортировка блоков пишет в
непересекающиеся отрезки `data_`. Слияние пишет в непересекающиеся отрезки
`merged_data`, читая `data_`. `data_.swap` выполняется после завершения
`parallel_for`.

## 6. Проверка корректности

Функциональный тест сравнивает результат с `std::ranges::sort`. Запуск текущего
`build_olesnitskiy/bin/ppc_func_tests` прошел для `seq/omp/stl/tbb`: 60 passed.
ALL отдельно прошел под `mpirun -np 2`: 15 passed.

## 7. Экспериментальная среда

- **Сборка:** `build_olesnitskiy`
- **Compiler:** `g++-14`
- **Flags:** `-O3 -DNDEBUG`, `std=gnu++23`
- **Размер входных данных:** `N=100000`
- **Baseline TaskRun:** `0.0058254364 s`
- **Baseline pipeline:** `0.0068995056 s`
- **Число повторов:** 5 по умолчанию

## 8. Результаты

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

## 9. Выводы

В измерениях TBB показал лучший результат среди потоковых backend-ов:
`0.0011774682 s` при 4 потоках, speedup `4.947`. Числа относятся к
`N=100000`; для других размеров нужны отдельные замеры.
