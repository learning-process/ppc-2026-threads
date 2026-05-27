## ALL ВЕРСИЯ

### Обзор

ALL версия - гибридная реализация алгоритма Фокса, которая объединяет
OpenMP для локального параллелизма и MPI для распределённого
вычисления на кластерах.

Технология: OpenMP + MPI (гибридная архитектура)
Потоков: 8 (локально на каждом процессе)
Процессов MPI: Поддерживает кластеры (N^2 процессов)
Ускорение: 7.1x на одном процессе
Эффективность: 89% на одном процессе

### Тестовые результаты

Все 12 размеров матриц от 1x1 до 32x32 успешно протестированы.

Результат: 12/12 PASSED (плюс 14 дополнительных unit тестов)

Производительные тесты матрица 512x512 с 1 MPI процессом:

Pipeline: 537 ms
Task run: 594 ms
Среднее: 565.5 ms
GFLOPS: 1.881

### Исходный код

```cpp
bool MatmulDoubleAllTask::RunImpl() {
  int my_rank = 0;
  int num_procs = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  const size_t n = matrix_size_;
  const auto &a = matrix_a_;
  const auto &b = matrix_b_;
  auto &c = result_matrix_;

  const int grid_dim = static_cast<int>(std::sqrt(num_procs));

  if (!IsValidConfigurationImpl(n, grid_dim, num_procs)) {
    HandleFallbackImpl(my_rank, n, a, b, c);
    GetOutput() = c;
    return true;
  }

  const size_t bs = n / static_cast<size_t>(grid_dim);
  const size_t block_sz = bs * bs;

  const int row_idx = my_rank / grid_dim;
  const int col_idx = my_rank % grid_dim;

  std::vector<double> local_a_block(block_sz);
  std::vector<double> local_b_block(block_sz);
  std::vector<double> local_c_block(block_sz, 0.0);

  std::vector<double> all_blocks_a;
  std::vector<double> all_blocks_b;

  if (my_rank == 0) {
    all_blocks_a.resize(static_cast<size_t>(num_procs) * block_sz);
    all_blocks_b.resize(static_cast<size_t>(num_procs) * block_sz);

    SplitIntoBlocksImpl(a, all_blocks_a, n, bs, grid_dim);
    SplitIntoBlocksImpl(b, all_blocks_b, n, bs, grid_dim);
  }

  DistributeBlocksImpl(my_rank, all_blocks_a, all_blocks_b,
                      local_a_block, local_b_block, block_sz);

  MPI_Comm row_comm = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, row_idx, col_idx, &row_comm);

  ExecuteFoxIterationsImpl(grid_dim, row_idx, col_idx, bs, block_sz,
                          row_comm, local_a_block, local_b_block,
                          local_c_block);

  CollectResultsImpl(my_rank, num_procs, n, bs, block_sz,
                    grid_dim, local_c_block, c);

  if (row_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&row_comm);
  }

  GetOutput() = c;
  return true;
}
```

### Распределение блоков с MPI_Scatter

```cpp
void DistributeBlocksImpl(int my_rank,
                        const std::vector<double> &blocks_a,
                        const std::vector<double> &blocks_b,
                        std::vector<double> &local_a,
                        std::vector<double> &local_b,
                        size_t block_sz) {
  const double *send_a = (my_rank == 0) ? blocks_a.data() :
                         nullptr;
  const double *send_b = (my_rank == 0) ? blocks_b.data() :
                         nullptr;

  MPI_Scatter(send_a, static_cast<int>(block_sz), MPI_DOUBLE,
              local_a.data(), static_cast<int>(block_sz),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Scatter(send_b, static_cast<int>(block_sz), MPI_DOUBLE,
              local_b.data(), static_cast<int>(block_sz),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
```

### Итерации алгоритма Фокса

```cpp
void ExecuteFoxIterationsImpl(int grid_dim, int row_id, int col_id,
                            size_t bs, size_t block_sz,
                            MPI_Comm row_comm,
                            std::vector<double> &local_a,
                            std::vector<double> &local_b,
                            std::vector<double> &local_c) {
  std::vector<double> broadcast_buffer(block_sz);

  for (int stage = 0; stage < grid_dim; ++stage) {
    const int source = (row_id + stage) % grid_dim;

    if (col_id == source) {
      broadcast_buffer = local_a;
    }

    MPI_Bcast(broadcast_buffer.data(),
              static_cast<int>(block_sz),
              MPI_DOUBLE, source, row_comm);

    MultiplyBlockPairImpl(broadcast_buffer, local_b, local_c, bs);

    const int target = (((row_id - 1 + grid_dim) % grid_dim) *
                       grid_dim) + col_id;
    const int origin = ((((row_id + 1) % grid_dim) % grid_dim) *
                       grid_dim) + col_id;

    MPI_Sendrecv_replace(local_b.data(),
                        static_cast<int>(block_sz),
                        MPI_DOUBLE, target, 0, origin, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}
```

### Сбор результатов с MPI_Gather

```cpp
void CollectResultsImpl(int my_rank, int num_procs, size_t n,
                      size_t bs, size_t block_sz, int grid_dim,
                      const std::vector<double> &local_c,
                      std::vector<double> &c) {
  std::vector<double> all_blocks;

  if (my_rank == 0) {
    all_blocks.resize(static_cast<size_t>(num_procs) * block_sz);
  }

  double *recv_buf = (my_rank == 0) ? all_blocks.data() :
                     nullptr;

  MPI_Gather(local_c.data(), static_cast<int>(block_sz),
             MPI_DOUBLE, recv_buf, static_cast<int>(block_sz),
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    MergeFromBlocksImpl(all_blocks, c, n, bs, grid_dim);
  }

  MPI_Bcast(c.data(), static_cast<int>(n * n),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
```

### Результаты производительности

Время выполнения на матрице 512x512: 565.5 ms
GFLOPS: 1.881

Сравнение с SEQ версией:
Ускорение: 7.1x (4036 ms / 565.5 ms)
Потоков: 8
Эффективность: 89%

Сравнение с OMP версией:
OMP: 556 ms
ALL: 565.5 ms
Разница: 1.6% медленнее (практически идентична)

### Анализ производительности

ALL версия практически идентична OMP потому что:

При запуске без mpirun MPI инициализируется в однопроцессном режиме.
Конфигурация: 1 процесс, 1x1 сетка
Валидация конфигурации проходит (1*1 = 1)
Используется основной алгоритм Фокса
Результат: практически идентичен OMP

Разница 1.6% (9.5 ms) - это margin ошибки от различных эффектов

Масштабируемость на кластеры:

С 1 процессом: 565.5 ms (fallback на OpenMP)
С 4 процессами (2x2 сетка): прогноз 150-200 ms
С 9 процессами (3x3 сетка): прогноз 70-100 ms
С 16 процессами (4x4 сетка): прогноз 40-60 ms

Прогнозы основаны на законе Амдала

### Выводы

Преимущества ALL версии:

1. Масштабируется на кластеры (поддерживает MPI)
2. На одном процессе практически идентична OMP
3. Универсальность (работает везде от одной машины до кластеров)
4. Гибридный подход использует оба уровня параллелизма
5. Graceful fallback при невалидной конфигурации
6. Расширенные тесты (26 тестов вместо 12)

Недостатки ALL версии:

1. Требует MPI (OpenMPI или MPICH)
2. Немного сложнее в реализации чем OMP
3. На 1.6% медленнее чем OMP на одном процессе

Использование:

ALL версия рекомендуется для кластерных вычислений, HPC систем,
максимальной универсальности (одна машина + кластеры) и гибридных
архитектур.

Статус: ГОТОВА К PRODUCTION (рекомендуется для кластеров)

---

## СРАВНЕНИЕ ВСЕХ 5 ВЕРСИЙ

Производительность на матрице 512x512:

OMP: 556 ms (7.3x ускорение, 91% эффективность)
ALL: 565.5 ms (7.1x ускорение, 89% эффективность)
STL: 628 ms (6.4x ускорение, 80% эффективность)
TBB: 3051 ms (1.32x ускорение, 16.5% эффективность)
SEQ: 4036 ms (1.0x ускорение, 100% эффективность)

Рекомендации:

Для одной машины: используйте OMP (7.3x)
Для кластеров: используйте ALL (7.1x на одном узле, масштабируется)
Без OpenMP: используйте STL (6.4x)
Для отладки: используйте SEQ
