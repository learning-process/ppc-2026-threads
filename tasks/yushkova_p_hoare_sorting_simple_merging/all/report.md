# Сортировка Хоара с простым слиянием - ALL

- **Студент:** Юшкова Полина Александровна, 3823Б1ПР2
- **Технология:** ALL (MPI + STL)
- **Вариант:** 13

## 1. Контекст

ALL-версия объединяет MPI между процессами и `std::thread` внутри каждого процесса. MPI используется для распределения
частей массива между rank-ами, а локальная STL-схема сортирует фрагмент внутри процесса.

## 2. Постановка задачи

- **Входные данные:** непустой `std::vector<int>`, доступный на корневом rank.
- **Выходные данные:** отсортированный по неубыванию `std::vector<int>`.
- **Baseline:** последовательная версия с временем `T_seq(task_run) = 0.0030497600 s`,
  `T_seq(pipeline) = 0.0083235000 s`.

## 3. Базовый алгоритм

Корневой rank хранит входной массив. Данные распределяются между rank-ами, каждый rank локально сортирует свой фрагмент
STL-схемой (блоки по 64 и простое слияние по уровням), затем фрагменты собираются на rank 0 и последовательно досливаются
в один общий массив. После этого итоговый вектор рассылается всем процессам.

## 4. Межпроцессная схема

`MPI_Comm_rank` и `MPI_Comm_size` определяют роль процесса и число rank-ов. `BuildDistribution` делит `total_size` почти
поровну: первые `remainder` rank-ов получают на один элемент больше. `MPI_Scatterv` отправляет локальные куски,
`MPI_Gatherv` собирает отсортированные куски на rank 0, затем `BroadcastVector` рассылает итог всем rank-ам через
`MPI_Bcast`.

Фрагмент распределения и обмена:

```cpp
std::vector<size_t> chunk_sizes(static_cast<size_t>(mpi_size));
std::vector<size_t> offsets(static_cast<size_t>(mpi_size));
BuildDistribution(total_size, mpi_size, chunk_sizes, offsets);

const std::vector<int> send_counts = MakeIntVector(chunk_sizes);
const std::vector<int> send_offsets = MakeIntVector(offsets);

std::vector<int> local_data(chunk_sizes[static_cast<size_t>(rank)]);
MPI_Scatterv(rank == 0 ? GetOutput().data() : nullptr, send_counts.data(),
             send_offsets.data(), MPI_INT, local_data.data(),
             send_counts[static_cast<size_t>(rank)], MPI_INT, 0,
             MPI_COMM_WORLD);

SortLocalStlParallel(local_data);

std::vector<int> gathered_data;
if (rank == 0) {
  gathered_data.resize(total_size);
}
MPI_Gatherv(local_data.data(), static_cast<int>(local_data.size()), MPI_INT,
            rank == 0 ? gathered_data.data() : nullptr, send_counts.data(),
            send_offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
```

## 5. Внутрипроцессная схема

Локальная сортировка выполняется функцией `SortLocalStlParallel`:

- разбиение на блоки по 64 элемента;
- сортировка каждого блока quicksort Хоара;
- уровни простого слияния по независимым парам диапазонов.

Задачи на каждом уровне выполняются через `RunInThreads` с распределением `task_index += thread_count`. Число потоков
выбирается по `std::thread::hardware_concurrency()` и числу задач (если `hardware_concurrency()==0`, используется
fallback `2`).

В коде явного `MPI_Barrier` нет; синхронизация возникает на коллективных операциях `MPI_Scatterv`, `MPI_Gatherv` и
`MPI_Bcast`.

## 6. Детали реализации

`ValidationImpl` проверяет непустой вход. `PreProcessingImpl` копирует вход в выходной буфер. `RunImpl` выполняет
MPI-распределение, локальную сортировку, сбор, финальное последовательное слияние на rank 0 и broadcast результата.
`PostProcessingImpl` проверяет, что выход непустой и отсортирован.

Файлы реализации: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.

## 7. Проверка корректности

ALL-backend зарегистрирован в общем функциональном тесте. Корректность подтверждается сравнением с эталоном
`std::ranges::sort` и дополнительными проверками `std::ranges::is_sorted`.

## 8. Экспериментальная среда

- **ОС:** Windows
- **MPI:** Microsoft MPI (`mpiexec`)
- **Compiler:** clang++ (MSVC toolchain), C++23
- **Размер входных данных:** `N=100000`
- **Диапазон значений:** `[-1000000, 1000000]`
- **Baseline TaskRun:** `0.0030497600 s`
- **Baseline pipeline:** `0.0083235000 s`
- **Число повторов:** 5

## 9. Результаты

- backend: all; ranks: 1; threads_per_rank: 12; total_workers: 12; time: 0.0159865600 s; speedup: 0.191;
  efficiency: 0.016; notes: `mpiexec -n 1`, local STL auto-threads, `TaskRun`.
- backend: all; ranks: 2; threads_per_rank: 12; total_workers: 24; time: 0.0155627800 s; speedup: 0.196;
  efficiency: 0.008; notes: `mpiexec -n 2`, local STL auto-threads, `TaskRun`.
- backend: all; ranks: 4; threads_per_rank: 12; total_workers: 48; time: 0.0125462200 s; speedup: 0.243;
  efficiency: 0.005; notes: `mpiexec -n 4`, local STL auto-threads, `TaskRun`.
- backend: all; ranks: 1; threads_per_rank: 12; total_workers: 12; time: 0.0173728400 s; speedup: 0.479;
  efficiency: 0.040; notes: `mpiexec -n 1`, local STL auto-threads, `pipeline`.
- backend: all; ranks: 4; threads_per_rank: 12; total_workers: 48; time: 0.0142562000 s; speedup: 0.584;
  efficiency: 0.012; notes: `mpiexec -n 4`, local STL auto-threads, `pipeline`.

## 10. Выводы

ALL-версия добавляет стоимость `Scatterv/Gatherv/Bcast` и финальное слияние на rank 0. На `N=100000` итоговое время
существенно зависит от коммуникаций и накладных расходов на создание локальных потоков, поэтому эффективность по
`total_workers` получается низкой, даже если локальная сортировка внутри ранга выполняется параллельно.
