# Умножение разреженных комплексных матриц (CCS) — ALL

- Student: Иманов Сабутай Ширзад оглы, группа 3823Б1ПР5
- Technology: ALL
- Variant: 7

## 1. Контекст

Технология **ALL** объединяет несколько средств курса в одном классе `SabutaySparseComplexCcsMultAll`:

1. **MPI** — разные процессы могут вызывать разные варианты умножения (SEQ, OMP, STL, TBB).
2. После умножения — короткие учебные примеры с OpenMP, `std::thread` и TBB.
3. В конце — `MPI_Barrier`, чтобы все процессы дождались друг друга.

В perf-тесте обычно один MPI-процесс (`PPC_NUM_PROC = 1`). Эталон SEQ: $T_{seq} = 0.00001390$ с.

## 2. Постановка задачи

Нужно вычислить произведение $C = A \cdot B$ для двух разреженных **комплексных** матриц в формате **CCS** (Compressed
Column Storage).

- **Вход:** `InType = std::tuple<CCS, CCS>` — матрицы $A$ и $B$.
- **Выход:** `OutType = CCS` — матрица $C$.
- **Ограничение на размеры:** $\text{cols}(A) = \text{rows}(B)$; иначе `ValidationImpl` возвращает ошибку.
- **Проверка структуры CCS:** длины `col_start`, `row_index`, `nz` согласованы; индексы строк в допустимом диапазоне.
- **Порог при слиянии:** элементы с модулем $\le 10^{-14}$ отбрасываются (`kDropMagnitude` в `all/src/ops_all.cpp`).
- **Критерий корректности в тестах:** совпадение с плотным умножением, допуск $10^{-12}$.

В технологии ALL основное умножение вызывается через `BuildProductMatrix`, а в `RunImpl` дополнительно демонстрируются
MPI, OpenMP, `std::thread` и TBB.

## 3. Базовый алгоритм

**Умножение (столбцовая схема CCS):**

1. для каждого столбца $j$ матрицы $B$ собрать вклад в столбец $C$ из столбцов $A$ (`BuildColumnFromRight`);
2. отсортировать пары (номер строки, значение) по строке;
3. объединить одинаковые строки и записать столбец в результирующую CCS.

**Выбор реализации умножения по MPI-рангу** (`rank % 4`):

| Остаток `rank % 4` | Backend |
| :----------------: | :-----: |
| 0 | SEQ (`SpmmAbc`) |
| 1 | OMP (`BuildProductMatrixOmp`) |
| 2 | STL (`SpmmAbc`) |
| 3 | TBB (`BuildProductMatrixTbb`) |

При `PPC_NUM_PROC = 1` всегда `rank = 0`, то есть умножение идёт как **SEQ**.

После умножения в `RunImpl` выполняются учебные фрагменты с потоками (не влияют на математику $C = A \cdot B$, но влияют
на время perf).

## 4. Схема работы

**Между процессами (MPI):**

- `MPI_Comm_rank` — номер процесса.
- `rank % 4`: 0 → SEQ, 1 → OMP, 2 → STL, 3 → TBB.
- При одном процессе всегда используется SEQ.

**Внутри одного процесса (после умножения):**

- на процессе 0 — короткий фрагмент OpenMP;
- затем создаются потоки `std::thread` (сразу `join`);
- затем `tbb::parallel_for` по счётчику.

### Листинг: выбор варианта и барьер

```287:325:tasks/sabutay_sparse_complex_ccs_mult_all/all/src/ops_all.cpp
bool SabutaySparseComplexCcsMultAll::RunImpl() {
  ...
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ppc::task::TypeOfTask backend = ppc::task::TypeOfTask::kSEQ;
  if (rank % 4 == 1) {
    backend = ppc::task::TypeOfTask::kOMP;
  } else if (rank % 4 == 2) {
    backend = ppc::task::TypeOfTask::kSTL;
  } else if (rank % 4 == 3) {
    backend = ppc::task::TypeOfTask::kTBB;
  }
  BuildProductMatrix(left, right, GetOutput(), backend);
  ...
  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}
```

## 5. Детали реализации

Основной код — `all/src/ops_all.cpp`, класс `SabutaySparseComplexCcsMultAll`.

**Конфигурация для perf-замеров в отчёте:**

| Параметр | Значение |
| -------- | -------- |
| MPI-процессы (`PPC_NUM_PROC`) | 1 |
| Потоков на процесс (`PPC_NUM_THREADS`) | 1, 2, 4, 8 |
| Умножение при 1 процессе | SEQ (`rank = 0`) |
| Эффективность $E$ | $S / p$, где $p$ = `PPC_NUM_THREADS` (как в остальных отчётах) |

После умножения — учебные фрагменты (OpenMP на `rank == 0`, цикл `std::thread`, `tbb::parallel_for`), затем
`MPI_Barrier`. Они **увеличивают** время perf по сравнению с чистым SEQ.

## 6. Проверка корректности

`tests/functional/main.cpp`, класс `SabutaySparseComplexCcsMultAll`.

**Результаты func-тестов:** 3/3 **PASSED**.

| `case_id` | Суффикс GTest | Результат |
| :-------: | ------------- | :-------: |
| 0 | `..._all_enabled_0` | PASSED |
| 1 | `..._all_enabled_1` | PASSED |
| 2 | `..._all_enabled_2` | PASSED |

```powershell
$env:PPC_NUM_PROC = "1"
mpiexec -n 1 build\bin\ppc_func_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*all_enabled*"
```

## 7. Экспериментальная среда

- **CPU:** Intel Core Ultra 5 125H, 14 ядер / 18 потоков.
- **RAM:** 32 ГБ, **ОС:** Windows 11, **сборка:** Release (MSVC).
- **Размер задачи (perf):** случайные матрицы $80 \times 90$ и $90 \times 70$.
- **MPI:** `PPC_NUM_PROC = 1`; **потоки:** `PPC_NUM_THREADS` = 1, 2, 4, 8.
- **Прогонов внутри замера:** 5.

```powershell
$env:PPC_NUM_PROC = "1"
$env:PPC_NUM_THREADS = "4"
build\bin\ppc_perf_tests.exe --gtest_filter="*sabutay_sparse_complex_ccs_mult_all*all_enabled*"
```

## 8. Результаты

Каркас perf запускает **`task_run`** и **`pipeline`** (1 MPI-процесс × $p$ потоков). $T_{seq} = 0.00001390$ с.

### `task_run`

| $p$ | Время (с) | $S$ | $E$ |
| :-: | :-------: | :-: | :-: |
| 1 | 0.00030686 | 0.05 | 5% |
| 2 | 0.00048390 | 0.03 | 1% |
| 4 | 0.00085534 | 0.02 | 0.5% |
| 8 | 0.00244152 | 0.01 | 0.1% |

_Комментарий:_ ALL **медленнее** чистого SEQ, потому что в `RunImpl` есть лишняя работа (демо-потоки и барьер MPI). Чем
больше `PPC_NUM_THREADS`, тем дольше идут эти дополнительные циклы.

### `pipeline`

| $p$ | Время (с) |
| :-: | :-------: |
| 1 | 0.00237558 |
| 2 | 0.00256790 |
| 4 | 0.00320222 |
| 8 | 0.00456848 |

## 9. Выводы

ALL — «сборная» версия для курса: MPI + несколько API. Для оценки скорости самого умножения смотрите отчёты SEQ, OMP и
TBB. ALL показывает, как связать разные технологии в одной программе.
