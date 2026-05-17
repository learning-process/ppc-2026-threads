# Умножение матриц: алгоритм Кэннона — ALL (MPI + OpenMP)

- Student: Зязева Светлана Александровна
- Technology: ALL (MPI + OpenMP)
- Variant: 1

## 1. Контекст

Гибридная версия реализует двухуровневый параллелизм: процессный уровень
(MPI) и внутрипроцессный уровень (OpenMP). Алгоритм Кэннона применяется
на уровне MPI-процессов — каждый процесс владеет одним блоком матрицы.
Внутри каждого процесса OpenMP дополнительно ускоряет вспомогательные
операции. Такая схема масштабируется как по числу узлов (MPI), так и по
числу ядер на каждом узле (OpenMP).

## 2. Постановка задачи

Аналогична SEQ: произведение квадратных матриц A и B размера N×N.
Baseline — `seq/report.md`. Конфигурация: `ranks x threads`, где число
MPI-рангов (`mpi_size_`) должно быть полным квадратом, а N должен
делиться на `sqrt(mpi_size_)`.

## 3. Базовый алгоритм

На уровне MPI реализован классический алгоритм Кэннона:

1. Ранг 0 распределяет блоки A и B по процессам (`DistributeBlocks`).
2. Каждый процесс выполняет начальное смещение блока A (по строке) и
   блока B (по столбцу) через `MPI_Sendrecv_replace`.
3. Выполняются `grid` шагов: умножение локальных блоков и циклический
   сдвиг A влево, B вверх через `MPI_Sendrecv_replace`.
4. Ранг 0 собирает результирующие блоки (`CollectResult`) и рассылает
   итоговую матрицу через `MPI_Bcast`.

При невозможности применить Кэннона ранг 0 выполняет
`RegularMultiplication` с OpenMP и рассылает результат.

## 4. Схема распараллеливания

**Процессный уровень (MPI):**

| Операция | Примитив |
| --- | --- |
| Рассылка размера матрицы | `MPI_Bcast` |
| Рассылка данных матриц | `MPI_Bcast` |
| Распределение блоков | `MPI_Send` / `MPI_Recv` |
| Начальное смещение A | `MPI_Sendrecv_replace` (по строке) |
| Начальное смещение B | `MPI_Sendrecv_replace` (по столбцу) |
| Cannon-шаги, сдвиг A | `MPI_Sendrecv_replace`, tag 20 |
| Cannon-шаги, сдвиг B | `MPI_Sendrecv_replace`, tag 21 |
| Сбор результата | `MPI_Recv` / `MPI_Send`, tag 30 |
| Рассылка результата | `MPI_Bcast` |

```cpp
// File: all/src/ops_all.cpp — Cannon-шаг
for (int step = 0; step < grid; ++step) {
  MultiplyBlocks(local_a, local_b, local_c, block_size);

  MPI_Sendrecv_replace(
      local_a.data(), block_elems, MPI_DOUBLE,
      left, 20, right, 20,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Sendrecv_replace(
      local_b.data(), block_elems, MPI_DOUBLE,
      up, 21, down, 21,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
```

`MPI_Sendrecv_replace` выполняет отправку и приём в один буфер,
исключая deadlock при кольцевом сдвиге. Все процессы вызывают функцию
одновременно — неявная синхронизация.

**Внутрипроцессный уровень (OpenMP):**

OpenMP применяется в:

- `RegularMultiplication` — `parallel for` по строкам (fallback);
- `AlignBlocks` — `parallel for collapse(2)` по блокам;
- `AssembleResult` — `parallel for collapse(2)` по блокам.

Везде используется `default(none)` с явным перечислением атрибутов.

**Валидация** выполняется только на ранге 0, результат рассылается:

```cpp
// File: all/src/ops_all.cpp — ValidationImpl
if (rank_ == 0) {
  valid = (sz > 0
           && m1.size() == sz * sz
           && m2.size() == sz * sz) ? 1 : 0;
}
MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
return valid != 0;
```

Это гарантирует согласованный результат валидации на всех рангах.

**Конфигурация `ranks x threads`:** реализация поддерживает произвольное
произведение. `PPC_NUM_PROC` (число MPI-процессов) и `PPC_NUM_THREADS`
(число OMP-потоков на процесс) задаются runner-ом курса.

## 5. Детали реализации

**Файлы:** `all/include/ops_all.hpp`, `all/src/ops_all.cpp`

`DistributeBlocks`: ранг 0 вычисляет координаты блока для каждого
процесса, копирует данные в буферы `tmp_a`/`tmp_b` и рассылает через
`MPI_Send`. Ранг 0 сохраняет свой блок напрямую.

```cpp
// File: all/src/ops_all.cpp — DistributeBlocks (упрощённо)
if (rank_ == 0) {
  for (int proc = 0; proc < mpi_size_; ++proc) {
    if (proc == 0) {
      local_a = tmp_a;
      local_b = tmp_b;
    } else {
      MPI_Send(tmp_a.data(), block_elems, MPI_DOUBLE,
               proc, 0, MPI_COMM_WORLD);
      MPI_Send(tmp_b.data(), block_elems, MPI_DOUBLE,
               proc, 1, MPI_COMM_WORLD);
    }
  }
} else {
  MPI_Recv(local_a.data(), block_elems, MPI_DOUBLE,
           0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(local_b.data(), block_elems, MPI_DOUBLE,
           0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
```

Замечание: последовательные `MPI_Send` на ранге 0 — узкое место при
большом числе процессов; альтернатива — `MPI_Scatterv`.

`CollectResult` симметрична: ранг 0 принимает блоки через `MPI_Recv` и
собирает результирующую матрицу.

## 6. Проверка корректности

ALL-версия сравнивалась с SEQ для всех функциональных тестов при числе
рангов 1, 4 (сетка 2×2) и при конфигурациях, активирующих fallback.
Синхронизированная валидация через `MPI_Bcast` гарантирует одинаковое
решение всех рангов. Расхождений не наблюдалось.

## 7. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | AMD Ryzen 7 8845HS w/ Radeon 780M Graphics |
| RAM | 32 ГБ |
| OS | Windows 11 Pro |
| Compiler | g++ / clang++ (Release) |
| CMake build type | Release |
| PPC_NUM_PROC | задаётся runner-ом |
| PPC_NUM_THREADS | задаётся runner-ом |

**Команды запуска:**

```bash
mpirun --allow-run-as-root -np 4 \
  ./build/bin/ppc_func_tests \
  --gtest_filter="*zyazeva_s_matrix_mult_cannon_alg_all*"

mpirun --allow-run-as-root -np 4 \
  ./build/bin/ppc_perf_tests \
  --gtest_filter="*zyazeva_s_matrix_mult_cannon_alg_all*"
```

## 8. Результаты

| Режим | Время ALL (с) | Время SEQ (с) | Ускорение |
| --- | --- | --- | --- |
| pipeline | 0.0720654500 | 0.7128969324 | **9.89x** |
| task_run | 0.0675401000 | 0.5624887756 | **8.33x** |

Гибридная версия показывает наибольшее ускорение среди всех реализаций.
Двухуровневый параллелизм позволяет задействовать несколько процессов и
все ядра внутри каждого процесса.

## 9. Выводы

MPI + OpenMP обеспечивают ускорение ~9.9x / ~8.3x. Ограничения: число
MPI-рангов должно быть полным квадратом, N — делиться на
`sqrt(mpi_size_)`; иначе fallback на ранге 0. Узкие места:
последовательная рассылка блоков через `MPI_Send` и финальный `MPI_Bcast`
всей матрицы, ограничивающий масштабируемость при больших N и большом
числе рангов.
