# Умножение плотных матриц. Элементы типа double. Блочная схема, алгоритм Кэннона. — ALL

- Student: Черемхин Андрей Александрович
- Technology: ALL
- Variant: 1

## 1. Контекст

ALL-версия реализует гибридную схему: распределение блоков между MPI-процессами и внутрипроцессное ускорение локальных операций с помощью OpenMP. 

## 2. Постановка задачи

На вход подаётся размер n и две матрицы A, B размера n x n в одномерном построчном виде. Требуется получить матрицу C = A * B. Валидация принимает только положительный n и два массива длины n * n.

Если n не совпадает с удобным блочным размером, реализация создаёт расширенные матрицы np x np и заполняет дополнительную область нулями. В финальный ответ копируется только исходная область n x n.

## 3. Базовый алгоритм

Матрица дополняется нулями до размера `padded_n`, который делится на размер виртуальной сетки. Каждый виртуальный узел сетки хранит блоки `A`, `B` и локальный блок результата `C`. На каждом шаге выполняется локальное умножение блоков, затем блоки `A` и `B` циклически сдвигаются по строкам и столбцам виртуальной сетки.

## 4. Межпроцессная схема

Реализация получает `world_rank` и `world_size` через MPI:

```cpp
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
```

Размер виртуальной сетки выбирается как минимальный `q`, для которого `q * q >= world_size`. Если виртуальных ячеек больше, чем реальных MPI-процессов, они распределяются циклически:

```cpp
int GetOwnerRank(int virtual_rank, int world_size) {
  return virtual_rank % world_size;
}
```

Rank 0 подготавливает padded-матрицы, извлекает из них начально сдвинутые блоки алгоритма Кэннона и отправляет их владельцам через `MPI_Send`. Остальные rank-и принимают свои блоки через `MPI_Recv`.

## 5. Внутрипроцессная схема

Внутри процесса используются OpenMP-циклы с `schedule(static)` для:

- копирования глобальной матрицы в padded-представление;
- извлечения и вставки локальных блоков;
- локального умножения блоков в `MulAddLocal`;
- копирования результата из padded-представления.

Число потоков задаётся через `PPC_NUM_THREADS` и применяется вызовом `omp_set_num_threads(requested_threads)`.

## 6. Детали реализации

Файлы реализации: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.

Основной цикл:

```cpp
for (int step = 0; step < q; ++step) {
  for (auto &cell : local_cells) {
    MulAddLocal(cell.a, cell.b, cell.c, block_n);
  }
  if (step + 1 < q) {
    ShiftBlocksCannon(local_cells, owner_by_rank, block_n, q, world_rank);
  }
}
```

Сдвиги выполняются в `ExchangePhase`. Если источник блока находится на том же MPI-процессе, выполняется локальное копирование. Если владелец другой, используются неблокирующие `MPI_Irecv` и `MPI_Isend`, затем `MPI_Waitall` для завершения обмена.

После вычисления `GatherResultBlocks` собирает блоки `C` на rank 0. Затем `MPI_Bcast` рассылает полную padded-матрицу результата всем rank-ам, чтобы каждый экземпляр задачи мог сформировать одинаковый `OutType`.

## 7. Проверка корректности

Корректность проверяется сравнением с oracle на тех же входах, что и у остальных backend-ов. Для гибридной версии важно запускать тесты при нескольких значениях `PPC_NUM_PROC`, потому что распределение виртуальных ячеек по владельцам меняется.

Команда:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Releas
cmake --build build -j --parallel
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# ranks = 2, threads per rank = 1, total workers = 2
PPC_NUM_PROC=2 PPC_NUM_THREADS=1 \
mpirun --oversubscribe -x PPC_NUM_PROC -x PPC_NUM_THREADS -n 2 \
  ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_all*'

# ranks = 2, threads per rank = 2, total workers = 4
PPC_NUM_PROC=2 PPC_NUM_THREADS=2 \
mpirun --oversubscribe -x PPC_NUM_PROC -x PPC_NUM_THREADS -n 2 \
  ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_all*'

# ranks = 4, threads per rank = 1, total workers = 4
PPC_NUM_PROC=4 PPC_NUM_THREADS=1 \
mpirun --oversubscribe -x PPC_NUM_PROC -x PPC_NUM_THREADS -n 4 \
  ./build/bin/ppc_perf_tests --gtest_filter='*cheremkhin_a_matr_mult_cannon_alg_all*'
```

Каждая команда `ppc_perf_tests` выводит две строки замеров: `task_run` и `pipeline`. Поэтому три конфигурации запуска соответствуют шести строкам таблицы результатов.

Дополнительно нужно проверять согласованность результата на всех rank-ах после `MPI_Bcast`.

## 8. Экспериментальная среда

- OS: Linux 6.6.114.1-microsoft-standard-WSL2 x86_64;
- CPU: AMD Ryzen 5 5600, 6 cores / 12 threads;
- Compiler: GCC 13.3.0 with MPI and OpenMP;
- Build type: `Release`;
- Process count: `PPC_NUM_PROC`;
- Threads per rank: `PPC_NUM_THREADS`;
- Normalization: `total_workers = PPC_NUM_PROC * PPC_NUM_THREADS`.

## 9. Результаты


| size | ranks | threads per rank | total workers | mode     | median time, s | speedup vs seq | efficiency |
| ---- | ----- | ---------------- | ------------- | -------- | -------------- | -------------- | ---------- |
| 640  | 2     | 1                | 2             | task     | 0.4014975420   | 6.0691         | 3.0346     |
| 640  | 2     | 2                | 4             | task     | 1.6332830144   | 1.4919         | 0.3730     |
| 640  | 4     | 1                | 4             | task     | 0.2357804242   | 10.3348        | 2.5837     |
| 640  | 2     | 1                | 2             | pipeline | 0.3924683468   | 6.2572         | 3.1286     |
| 640  | 2     | 2                | 4             | pipeline | 1.6682519290   | 1.4720         | 0.3680     |
| 640  | 4     | 1                | 4             | pipeline | 0.2316247760   | 10.6022        | 2.6506     |


Ускорение рассчитано как `T_seq / T_all`, а эффективность — как `speedup / total_workers`. Использованный SEQ baseline: `2.4367350032` для `task` и `2.4557342674` для `pipeline`. При интерпретации результатов нужно отделять вычислительный выигрыш от стоимости MPI-обменов. На размере `640 x 640` сдвиги блоков и финальный broadcast могут быть заметной частью времени.

## 10. Выводы

Гибридная версия показывает полный вариант алгоритма Кэннона с распределённой сеткой блоков. Она методически наиболее близка к исходной идее алгоритма, но её эффективность зависит от размера матрицы, числа rank-ов, числа потоков внутри rank-а и цены коммуникации.