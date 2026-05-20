# Сортировка Хоара с простым слиянием — ALL

- **Студент:** Олесницкий Владимир Тарасович, 3823Б1ПР2
- **Технология:** ALL (MPI + STL)
- **Вариант:** 13

## 1. Контекст

ALL-версия объединяет MPI между процессами и `std::thread` внутри каждого
процесса. MPI используется для распределения частей массива между rank-ами, а
локальная STL-схема сортирует фрагмент внутри процесса.

## 2. Постановка задачи

- **Входные данные:** непустой `std::vector<int>`, доступный на корневом rank.
- **Выходные данные:** отсортированный по неубыванию `std::vector<int>`.
- **Baseline:** последовательная версия с временем `T_seq = 0.0058254364 s`.

## 3. Базовый алгоритм

Корневой rank хранит входной массив. Данные распределяются между rank-ами,
каждый rank локально сортирует свой фрагмент STL-схемой, затем фрагменты
собираются на rank 0 и последовательно досливаются.

## 4. Межпроцессная схема

`MPI_Comm_rank` и `MPI_Comm_size` определяют роль процесса и число rank-ов.
`BuildDistribution` делит `total_size` почти поровну: первые `remainder`
rank-ы получают на один элемент больше. `MPI_Scatterv` отправляет локальные
куски, `MPI_Gatherv` собирает отсортированные куски на rank 0, затем
`BroadcastVector` рассылает итог всем rank-ам через `MPI_Bcast`.

Фрагмент распределения и обмена:

```cpp
std::vector<size_t> chunk_sizes(static_cast<size_t>(mpi_size));
std::vector<size_t> offsets(static_cast<size_t>(mpi_size));
BuildDistribution(total_size, mpi_size, chunk_sizes, offsets);

std::vector<int> send_counts = MakeIntVector(chunk_sizes);
std::vector<int> send_offsets = MakeIntVector(offsets);
std::vector<int> local_data(chunk_sizes[static_cast<size_t>(mpi_rank)]);

MPI_Scatterv(data_.data(), send_counts.data(), send_offsets.data(), MPI_INT,
             local_data.data(), send_counts[mpi_rank], MPI_INT, 0,
             MPI_COMM_WORLD);

SortLocalStlParallel(local_data);

std::vector<int> gathered_data;
if (mpi_rank == 0) {
  gathered_data.resize(total_size);
}
```

## 5. Внутрипроцессная схема

Локальная сортировка повторяет STL-версию: `RunInThreads`, блоки по 64, затем
уровни слияния. Количество потоков внутри rank-а выбирается по
`std::thread::hardware_concurrency()` и числу локальных задач.

Для ALL важно указывать `ranks`, `threads_per_rank` и
`total_workers = ranks * threads_per_rank`. В коде явного `MPI_Barrier` нет;
обязательная синхронизация возникает на коллективных операциях `MPI_Scatterv`,
`MPI_Gatherv` и `MPI_Bcast`. Дополнительный `MPI_Barrier` добавляет тестовый
listener раннера после каждого теста.

## 6. Детали реализации

`ValidationImpl` проверяет непустой вход. `PreProcessingImpl` копирует вход в
`data_`. `RunImpl` выполняет MPI-распределение, локальную сортировку, сбор и
broadcast. `PostProcessingImpl` проверяет `std::ranges::is_sorted` и записывает
выход.

## 7. Проверка корректности

ALL-backend зарегистрирован в общем функциональном тесте. Запуск под
`mpirun -np 2` с фильтром
`*olesnitskiy_v_hoare_sort_simple_merge_all_enabled*` прошел 15/15 ALL-тестов.
Остальные backend-ы (`seq/omp/stl/tbb`) прошли 60 тестов вне `mpirun`.

## 8. Экспериментальная среда

- **Сборка:** `build_olesnitskiy`
- **Compiler:** `g++-14`
- **Flags:** `-O3 -DNDEBUG`, `std=gnu++23`
- **Размер входных данных:** `N=100000`
- **Baseline TaskRun:** `0.0058254364 s`
- **Baseline pipeline:** `0.0068995056 s`
- **Число повторов:** 5 по умолчанию

## 9. Результаты

- backend: all; ranks: 1; threads_per_rank: 12; total_workers: 12; time:
  0.0126647332 s; speedup: 0.460; efficiency: 0.038; notes: `mpirun -np 1`,
  `PPC_NUM_THREADS=1`, local STL auto-threads.
- backend: all; ranks: 2; threads_per_rank: 12; total_workers: 24; time:
  0.0088037862 s; speedup: 0.662; efficiency: 0.028; notes: `mpirun -np 2`,
  local STL auto-threads.
- backend: all; ranks: 4; threads_per_rank: 12; total_workers: 48; time:
  0.0043261284 s; speedup: 1.347; efficiency: 0.028; notes: `mpirun -np 4`,
  local STL auto-threads.
- backend: all; ranks: 4; threads_per_rank: 12; total_workers: 48; time:
  0.0516886980 s; speedup: 0.133; efficiency: 0.003; notes: `pipeline`,
  `mpirun -np 4`, local STL auto-threads.

## 10. Выводы

ALL-версия добавляет стоимость `Scatterv/Gatherv/Bcast` и финальное слияние на
rank 0. На 4 rank-ах получено `0.0043261284 s`, speedup `1.347`; efficiency
`0.028`, потому что каждый rank дополнительно запускает до 12 STL-потоков, а
коммуникации остаются обязательными.
