# Отчет ALL: MPI + STL для сортировки Хоара с простым слиянием

## Контекст и базовый алгоритм

ALL-версия объединяет MPI между процессами и `std::thread` внутри процесса. Корневой rank хранит входной массив, данные
распределяются между rank-ами, каждый rank локально сортирует свой фрагмент STL-схемой, затем фрагменты собираются на
rank 0 и последовательно досливаются ([`all/src/ops_all.cpp`](src/ops_all.cpp#L227)).

## Межпроцессная схема

`MPI_Comm_rank` и `MPI_Comm_size` определяют роль процесса и число rank-ов ([`all/src/ops_all.cpp`](src/ops_all.cpp#L228)).
`BuildDistribution` делит `total_size` почти поровну: первые `remainder` rank-ов получают на один элемент больше
([`all/src/ops_all.cpp`](src/ops_all.cpp#L64)). `MPI_Scatterv` отправляет локальные куски, `MPI_Gatherv` собирает
отсортированные куски на rank 0, затем `BroadcastVector` рассылает итог всем rank-ам через `MPI_Bcast`
([`all/src/ops_all.cpp`](src/ops_all.cpp#L246), [`all/src/ops_all.cpp`](src/ops_all.cpp#L256),
[`all/src/ops_all.cpp`](src/ops_all.cpp#L264)).

Фрагмент, [`all/src/ops_all.cpp`](src/ops_all.cpp#L238): распределение, scatter/gather и broadcast.

```cpp
std::vector<size_t> chunk_sizes(static_cast<size_t>(mpi_size));
std::vector<size_t> offsets(static_cast<size_t>(mpi_size));
BuildDistribution(total_size, mpi_size, chunk_sizes, offsets);

std::vector<int> send_counts = MakeIntVector(chunk_sizes);
std::vector<int> send_offsets = MakeIntVector(offsets);
std::vector<int> local_data(chunk_sizes[static_cast<size_t>(mpi_rank)]);

MPI_Scatterv(data_.data(), send_counts.data(), send_offsets.data(), MPI_INT, local_data.data(), send_counts[mpi_rank],
             MPI_INT, 0, MPI_COMM_WORLD);

SortLocalStlParallel(local_data);

std::vector<int> gathered_data;
if (mpi_rank == 0) {
  gathered_data.resize(total_size);
}
```

## Внутрипроцессная схема

Локальная сортировка повторяет STL-версию: `RunInThreads`, блоки по 64, затем уровни слияния
([`all/src/ops_all.cpp`](src/ops_all.cpp#L167)). Количество потоков внутри rank-а выбирается по
`std::thread::hardware_concurrency()` и числу локальных задач ([`all/src/ops_all.cpp`](src/ops_all.cpp#L20)).

## Конфигурация ranks × threads и цена синхронизации

Для ALL нужно указывать `ranks`, `threads_per_rank`, `total_workers = ranks * threads_per_rank`. В коде явного
`MPI_Barrier` нет; обязательная цена синхронизации возникает на коллективных операциях `MPI_Scatterv`, `MPI_Gatherv` и
`MPI_Bcast`. Дополнительный `MPI_Barrier` добавляет тестовый listener раннера после каждого теста
([`modules/runners/src/runners.cpp`](../../../modules/runners/src/runners.cpp#L23)).

## Детали pipeline и корректность

`ValidationImpl` проверяет непустой вход ([`all/src/ops_all.cpp`](src/ops_all.cpp#L217)). `PreProcessingImpl` копирует
вход в `data_` ([`all/src/ops_all.cpp`](src/ops_all.cpp#L221)). `RunImpl` выполняет MPI-распределение, локальную
сортировку, сбор и broadcast ([`all/src/ops_all.cpp`](src/ops_all.cpp#L227)). `PostProcessingImpl` проверяет
`std::ranges::is_sorted` и записывает выход ([`all/src/ops_all.cpp`](src/ops_all.cpp#L268)). В исходном функциональном
тесте ALL-backend зарегистрирован ([`tests/functional/main.cpp`](../tests/functional/main.cpp#L88)).

## Результаты

Baseline: `seq` `TaskRun = 0.0058254364 s`; для pipeline baseline `0.0068995056 s`. Performance-вход `N=100000`
([`tests/performance/main.cpp`](../tests/performance/main.cpp#L20)). Framework выполняет 5 повторов по умолчанию
([`modules/performance/include/performance.hpp`](../../../modules/performance/include/performance.hpp#L21)).

| backend | ranks | threads_per_rank | total_workers | time | speedup | efficiency | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| all | 1 | 12 | 12 | 0.0126647332 s | 0.460 | 0.038 | `mpirun -np 1`, `PPC_NUM_THREADS=1`, local STL auto-threads |
| all | 2 | 12 | 24 | 0.0088037862 s | 0.662 | 0.028 | `mpirun -np 2`, local STL auto-threads |
| all | 4 | 12 | 48 | 0.0043261284 s | 1.347 | 0.028 | `mpirun -np 4`, local STL auto-threads |
| all | 4 | 12 | 48 | 0.0516886980 s | 0.133 | 0.003 | `pipeline`, `mpirun -np 4`, local STL auto-threads |

## Выводы

ALL-версия добавляет стоимость `Scatterv/Gatherv/Bcast` и финальное слияние на rank 0
([`all/src/ops_all.cpp`](src/ops_all.cpp#L204)). На 4 rank-ах получено `0.0043261284 s`, speedup `1.347`; efficiency
`0.028`, потому что каждый rank дополнительно запускает до 12 STL-потоков, а коммуникации остаются обязательными.
