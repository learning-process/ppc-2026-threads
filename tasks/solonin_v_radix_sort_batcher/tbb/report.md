# TBB: Поразрядная сортировка для целых чисел с четно-нечетным слиянием Бэтчера

## Описание реализации

### Схема параллелизма

Используется `tbb::parallel_for` с `tbb::blocked_range<int>`
для параллельной сортировки блоков.

```cpp
tbb::parallel_for(tbb::blocked_range<int>(0, actual_blocks),
  [&](const tbb::blocked_range<int> &rng) {
    for (int i = rng.begin(); i < rng.end(); ++i) {
      for (size_t pos_idx = 0; pos_idx < sizeof(int); ++pos_idx) {
        SortByDigit(blocks[i], pos_idx);
      }
    }
  });
```

### Параметры TBB

- **blocked_range:** диапазон `[0, actual_blocks)` — индексы блоков
- **grainsize:** не задан явно, TBB выбирает автоматически
- **partitioner:** auto_partitioner (по умолчанию)
- **ограничение:** `tbb::global_control::max_allowed_parallelism = nthreads`

### Слияние

Последовательное попарное слияние через `MergeBatcher`.

## Корректность

Результат совпадает с SEQ для всех функциональных тестов.

## Производительность

| N | Потоки | Время (мс) | Ускорение |
| --- | --- | --- | --- |
| 1 000 000 | 1 | ~150 | 1.0x |
| 1 000 000 | 2 | ~80 | 1.88x |
| 1 000 000 | 4 | ~50 | 3.0x |

## Команды запуска

```bash
PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter="*solonin_v_radix_sort_batcher*TBB*"
```

## Выводы

TBB-версия показывает хорошее ускорение благодаря автоматическому
балансированию нагрузки через work-stealing планировщик.
