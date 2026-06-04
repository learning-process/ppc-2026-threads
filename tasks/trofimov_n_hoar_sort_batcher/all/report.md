# Hoare Sort + Batcher Merge — ALL (MPI + OMP + TBB + STL)

- Студент: Трофимов Никита Сергеевич
- Технология: ALL
- Вариант: 14

## Межпроцессная схема (MPI)

1. Rank 0 строит `counts/displs`.
2. `MPI_Scatterv` распределяет чанки.
3. Локальная сортировка в каждом процессе.
4. `MPI_Gatherv` собирает отсортированные части на rank 0.
5. Rank 0 сливает чанки и выполняет `MPI_Bcast` итогового массива.

## Внутрипроцессная схема

Локальные участки используют OMP/TBB/STL-компоненты для ускорения сортировки.

## Реализация

- Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`
- Этапы `RunImpl`: `BuildDistribution`, `ScatterData`, `SortLocalChunkHybrid`,
  `GatherChunks`, `MergeSortedChunksOnRoot`, `BroadcastResult`.

## Корректность

Проверяется глобальная отсортированность и согласованность результата между рангами.

## Производительность (из логов)

- task run: `0.002739 s`, speedup `5.70`
- pipeline: `0.010874 s`, speedup `3.87`

## Вывод

ALL объединяет все требуемые технологии и даёт высокий выигрыш,
но чувствителен к накладным расходам MPI.
