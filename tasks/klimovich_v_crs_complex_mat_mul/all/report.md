# Умножение разреженных матриц комплексного типа в формате CRS — ALL (MPI + OMP)

- Student: Климович Виктор Олегович, group 3823Б1ПР4
- Technology: MPI + OpenMP
- Variant: 6

## 1. Контекст

ALL-версия — гибридная: на верхнем уровне MPI делит строки матрицы `A`
между процессами (rank-ами), внутри каждого rank-а параллельный обход
строк выполняется через OpenMP. Это самый методически нагруженный отчёт,
потому что здесь нужно описать одновременно межпроцессную и
внутрипроцессную схемы, явно указывать коммуникационные паттерны и
обосновывать выбранную нормировку эффективности `ranks × threads`.

## 2. Постановка задачи

См. [seq/report.md](../seq/report.md). Дополнительное требование:
результат строится только на rank 0; на остальных rank-ах `GetOutput()`
после `RunImpl` остаётся в исходном default-состоянии.

## 3. Базовый алгоритм

Тот же Густавсон row-wise со SPA внутри каждого процесса. Новое — два
уровня распараллеливания: строки `[0, n_rows)` делятся между `world_size`
процессами, затем строки локального диапазона делятся между OMP-потоками
внутри процесса.

## 4. Межпроцессная схема

- Только rank 0 получает входные данные через конструктор
  (`if (rank == 0) GetInput() = in;`). Остальные ранги начинают с пустых
  `lhs`, `rhs`.
- На входе в `RunImpl` rank 0 рассылает обе матрицы через
  `BroadcastOperand`: сначала метаданные `(n_rows, n_cols, nnz)` массивом
  из 3 `int`, затем `row_offsets`, `col_indices` и `data`. Комплексные
  значения передаются как пары `double` (`nnz * 2` элементов), что
  корректно ровно потому, что layout `std::complex<double>` совпадает с
  двумя последовательными `double` по стандарту.
- Каждый rank вычисляет свой диапазон строк через `RowRange`:
  `base = n_rows / world_size`, `extra = n_rows % world_size`, первые
  `extra` rank-ов получают `base + 1` строк, остальные — `base`. Это та
  же формула, что в STL-версии, но на уровне процессов.
- После локального счёта каждый rank возвращает три потока данных: число
  ненулей на строку (`local_nnz_per_row`), столбцы (`local_cols`) и
  комплексные значения (`local_vals`). Сборка идёт через `MPI_Gatherv`:
  - `local_nnz_per_row` собирается с counts/displs, выраженными в
    **строках**;
  - `local_cols` собирается с counts/displs в **штуках элементов**;
  - `local_vals` использует те же counts/displs, **умноженные на 2**
    (потому что значения отправляются как `MPI_DOUBLE`).
- На rank 0 итоговый `row_offsets` строится последовательно префиксной
  суммой по `global_nnz_per_row`, а `col_indices`/`data` присваиваются из
  собранных векторов через `std::move`.

Используемые MPI-вызовы: `MPI_Comm_rank`, `MPI_Comm_size`, `MPI_Bcast`
(для рассылки операндов), `MPI_Gather` (для одиночных payload-counts),
`MPI_Gatherv` (для трёх payload-векторов). Барьеры явно не вызываются —
коммуникационные операции сами по себе синхронизируют участников.

## 5. Внутрипроцессная схема

Внутри `ComputeLocalRows` каждый rank запускает OMP-параллель ровно по
своему диапазону `[row_begin, row_end)`. Схема та же, что в OMP-версии:
приватные `spa`, `touched_by_row`, `touched_cols`,
`schedule(dynamic, 16)`. Это сознательное переиспользование рабочей
OMP-схемы, чтобы один и тот же параллельный код использовался и в чистом
OMP-backend-е, и внутри ALL.

## 6. Детали реализации

Файлы: [all/include/ops_all.hpp](include/ops_all.hpp),
[all/src/ops_all.cpp](src/ops_all.cpp).

Бродкаст операнда:

```cpp
// File: all/src/ops_all.cpp
void KlimovichVCrsComplexMatMulAll::BroadcastOperand(CrsMatrix &m, int root) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 3> meta{0, 0, 0};
  if (rank == root) {
    meta[0] = m.n_rows;
    meta[1] = m.n_cols;
    meta[2] = static_cast<int>(m.data.size());
  }
  MPI_Bcast(meta.data(), 3, MPI_INT, root, MPI_COMM_WORLD);

  if (rank != root) {
    m.n_rows = meta[0];
    m.n_cols = meta[1];
    m.row_offsets.assign(static_cast<std::size_t>(meta[0]) + 1, 0);
    m.col_indices.assign(static_cast<std::size_t>(meta[2]), 0);
    m.data.assign(static_cast<std::size_t>(meta[2]), Cplx(0.0, 0.0));
  }

  if (meta[0] > 0) {
    MPI_Bcast(m.row_offsets.data(), meta[0] + 1, MPI_INT, root, MPI_COMM_WORLD);
  }
  if (meta[2] > 0) {
    MPI_Bcast(m.col_indices.data(), meta[2], MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(m.data.data(), meta[2] * 2, MPI_DOUBLE, root, MPI_COMM_WORLD);
  }
}
```

Метаданные шлются всегда (`MPI_Bcast` с 3 `int`); тяжёлые массивы —
только если они непусты. На non-root ранке буферы предварительно
выделяются под известный из меты размер.

Параллельный счёт внутри rank-а:

```cpp
// File: all/src/ops_all.cpp
#pragma omp parallel default(none) shared(lhs, rhs, stages, row_begin, row_end)
{
  std::vector<Cplx> spa(static_cast<std::size_t>(rhs.n_cols));
  std::vector<int> touched_by_row(static_cast<std::size_t>(rhs.n_cols), -1);
  std::vector<int> touched_cols;
  touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));

#pragma omp for schedule(dynamic, 16)
  for (int i = row_begin; i < row_end; ++i) {
    GustavsonRow(lhs, rhs, i, spa, touched_by_row, touched_cols, stages[i - row_begin]);
  }
}
```

Та же `value-initialization` (`std::vector<Cplx> spa(n)` без второго
аргумента) применена и здесь — иначе MSVC выдаёт C3052 на
`Cplx(0.0, 0.0)` под `default(none)`.

Чтобы снизить cognitive complexity функции `RunImpl` (требование
clang-tidy), часть MPI-обвязки вынесена в свободные функции
`FillRowsPerProc`, `GatherPayloadCountsAndDispls`, `BuildPayloadCountsD`
— они выполняют только подготовку counts/displs для `MPI_Gatherv`.

## 7. Проверка корректности

Functional-тест [tests/functional/main.cpp](../tests/functional/main.cpp)
прогоняет ALL-версию через
`AddFuncTask<KlimovichVCrsComplexMatMulAll, InType>` под `mpirun` на тех
же 10 наборах размеров. Проверка результата выполняется только на rank 0
(`CheckTestOutputData` делает `rank != 0 → return true`). Все наборы
проходят с допуском `1e-9` для конфигураций `mpirun -n N` с
`N ∈ {1, 2, 4}`.

Согласованность результата между MPI-конфигурациями: при разных `N` итог
на rank 0 побайтно совпадает со SEQ-результатом, потому что разбиение
строк не меняет порядок суммирования внутри одной строки.

## 8. Экспериментальная среда

- CPU: Intel Core i7-11800H @ 2.30 GHz (8C/16T)
- RAM: 16 GB DDR4-3200
- OS: Windows 11 Pro 22H2 (build 22631)
- Compiler: MSVC 19.44.35211 (Visual Studio 2022)
- MPI: Microsoft MPI v10.1.12498.18
- Build type: Release
- `PPC_NUM_PROC` ∈ {1, 2, 4}; `PPC_NUM_THREADS` ∈ {1, 2, 4} (комбинации
  `ranks × threads`).

Команды:

```bash
set PPC_NUM_PROC=2
set PPC_NUM_THREADS=2
scripts/run_tests.py --running-type=processes --counts 2 4
scripts/run_tests.py --running-type=performance
```

## 9. Результаты

Размер задачи `1500×1500`, `T_seq(task) = 0.152 s`,
`T_seq(pipeline) = 1.498 s`. Медианы по 10 повторам.

Эффективность считается как
`E = T_seq / (T(ranks, threads) * total_workers)`, где
`total_workers = ranks * threads/rank`. Это та же нормировка, что
используется в OMP/TBB/STL по числу потоков — выбрана единая база, чтобы
сравнение технологий в корневом отчёте было методически корректным.

| ranks | threads/rank | total workers | mode     | median time, s | speedup vs seq | efficiency, % |
| ----: | -----------: | ------------: | -------- | -------------: | -------------: | ------------: |
|     1 |            1 |             1 | task     |          0.165 |           0.92 |            92 |
|     1 |            4 |             4 | task     |          0.052 |           2.92 |            73 |
|     2 |            2 |             4 | task     |          0.063 |           2.41 |            60 |
|     4 |            1 |             4 | task     |          0.085 |           1.79 |            45 |
|     4 |            2 |             8 | task     |          0.058 |           2.62 |            33 |
|     1 |            1 |             1 | pipeline |          1.610 |           0.93 |            93 |
|     1 |            4 |             4 | pipeline |          0.452 |           3.31 |            83 |
|     2 |            2 |             4 | pipeline |          0.594 |           2.52 |            63 |
|     4 |            1 |             4 | pipeline |          0.812 |           1.84 |            46 |
|     4 |            2 |             8 | pipeline |          0.540 |           2.77 |            35 |

Наблюдения:

- Лучшая конфигурация на одной машине — `ranks=1, threads=4`: фактически
  это чистый OMP-режим без MPI-overhead, и результат близок к чистому
  OMP-backend-у.
- Конфигурации с `ranks > 1` платят коммуникационным overhead:
  `BroadcastOperand` для обоих операндов выполняется каждый запуск,
  `MPI_Gatherv` собирает результат на rank 0.
- Чем больше rank-ов, тем больше cost broadcast-ов (отправка одних и тех
  же данных N-1 раз) и сильнее последовательное узкое место в Gatherv на
  rank 0 → эффективность падает с 60% (2×2) до 33% (4×2).

Цена коммуникации, цифры порядка величины: оба операнда — это
`≈ 2 * 1500 * (sizeof(int) + sizeof(double) * 2) = ~72 KB` на rank.
Это мизерно по сравнению с временем вычисления (десятки миллисекунд) на
shared-memory MPI; на сетевом MPI этот же бродкаст занял бы заметную
долю времени.

## 10. Выводы

Гибридная схема оправдана, когда на машине доступно больше суммарных
worker-ов, чем поддерживает один OMP-runtime эффективно (например,
NUMA-узлы или distributed-memory кластер). На одной машине с одним
сокетом `ranks=1, threads=N` обычно выигрывает у `ranks=N, threads=1` за
счёт отсутствия MPI-overhead на бродкасты и Gatherv — это и подтверждают
замеры (5.2× у OMP-8 против 2.6× у ALL-4×2). Главное методическое
решение этого backend-а — переиспользовать чистую OMP-схему внутри
каждого rank-а, что упрощает доказательство корректности и сокращает
количество новых синхронизационных точек до одной (`MPI_Gatherv` в
конце).
