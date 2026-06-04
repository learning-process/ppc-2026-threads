# Умножение разреженных матриц (CCS, double) — ALL (MPI + OpenMP)

- Student: Дилшодов Адхам Умидович, группа 3823Б1ПР4
- Technology: ALL — гибридная (MPI + OpenMP)
- Variant: 5

## 1. Контекст

Самая сложная из реализаций: вычисление одновременно распределяется между процессами (MPI)
и потоками внутри каждого процесса (OpenMP). Цель — продемонстрировать иерархический
параллелизм, в котором рабочая единица определяется парой `(ranks × threads)`,
а не одним числом потоков.

## 2. Постановка задачи

Без изменений: `C = A · B`, оба операнда в CCS, инварианты структуры и численная точность
те же, что в [seq/report.md](../seq/report.md). Особенность для гибридной версии: результат
должен быть доступен на **каждом** MPI-ранге после `RunImpl`, поскольку инфраструктура курса
вызывает `CheckTestOutputData(task_->GetOutput())` на всех рангах.

## 3. Базовый алгоритм

Тот же sparse-accumulator паттерн, но обёрнутый в два уровня декомпозиции:

1. **Между процессами:** диапазон столбцов `[0, cols_B)` равномерно разбивается на
   `world_size` непересекающихся «полос» (slab-ов);
2. **Внутри процесса:** полоса обрабатывается OpenMP-параллельным циклом по столбцам,
   каждый поток держит свои `acc`/`marker`/`used_rows`.

После локального вычисления выполняется `MPI_Allgatherv` для `row_indices` и `values`,
плюс отдельный `Allgatherv` для сдвинутых локальных `col_ptrs`.
Результат имеется на каждом ранге.

## 4. Межпроцессная схема

- **Коммуникатор:** `MPI_COMM_WORLD`. Каждый ранг получает `rank` и `size`
  через `MPI_Comm_rank/Comm_size`.
- **Роли рангов:** симметричные — нет управляющего процесса. Все ранги выполняют одно
  и то же распределение, считают свою полосу и участвуют в коллективных операциях.
- **Распределение данных:** в текущей реализации **обе матрицы** (`lhs`, `rhs`) уже находятся
  на каждом ранге, потому что инфраструктура курса вызывает `SetUp()` независимо на каждом
  MPI-ранке и `GetTestInputData()` возвращает полный вход. Поэтому отдельный `MPI_Bcast`
  входов не нужен. В сценарии «реальной» production-нагрузки потребовался бы
  `Bcast(A)` + `Scatter(B by columns)`.
- **Точки коллективной синхронизации:**
  1. `MPI_Allgather(&local_nnz)` — собрать число ненулей с каждого ранга;
  2. `MPI_Allgatherv(row_indices)` — объединить индексы строк;
  3. `MPI_Allgatherv(values)` — объединить значения;
  4. `MPI_Allgatherv(local_inner)` — объединить сдвинутые `col_ptrs`
     (см. ниже про сдвиг на `nnz_displs[rank]`).

Сдвиг `col_ptrs` — нетривиальный шаг. Локальные `col_ptrs` в каждом ранге начинаются с нуля
и индексируют **локальный** буфер ненулей. Чтобы превратить их в глобальные индексы,
перед `Allgatherv` каждый элемент `local_col_ptrs[c]` увеличивается на `nnz_displs[rank]` —
суммарное число ненулей всех предыдущих рангов:

```cpp
// File: all/src/ops_all.cpp
std::vector<int> local_inner(local_cols);
for (int col_idx = 0; col_idx < local_cols; ++col_idx) {
  local_inner[col_idx] = local_col_ptrs[col_idx] + nnz_displs[rank];
}
MPI_Allgatherv(local_inner.data(), local_cols, MPI_INT,
               out.col_ptrs.data(), cols_per_proc.data(), col_displs.data(),
               MPI_INT, MPI_COMM_WORLD);
out.col_ptrs[out.cols_count] = total_nnz; // последний элемент
```

После этой операции `out.col_ptrs` содержит глобальные смещения, согласованные
с конкатенированными `row_indices`/`values`.

## 5. Внутрипроцессная схема

Внутри ранга работа разворачивается OpenMP-ом по тому же шаблону, что и в чистом OMP-backend:

```cpp
// File: all/src/ops_all.cpp
#pragma omp parallel default(none) shared(lhs, rhs, column_results, j_start, j_end)
{
  std::vector<double> acc(lhs.rows_count, 0.0);
  std::vector<int>    marker(lhs.rows_count, -1);
  std::vector<int>    used_rows;

#pragma omp for schedule(dynamic)
  for (int j = j_start; j < j_end; ++j) {
    AccumulateColumnProduct(lhs, rhs, j, acc, marker, used_rows, column_results[j - j_start]);
  }
}
```

- Каждый OMP-поток держит приватные `acc`, `marker`, `used_rows`;
- `schedule(dynamic)` — устойчивость к неравномерной плотности колонок;
- `column_results` индексируется **локальным** номером колонки (`j - j_start`), поэтому потоки
  разных рангов и потоки внутри одного ранга не пересекаются по записи.

После выхода из `parallel`-региона делается последовательное построение `local_col_ptrs`,
`local_row_indices`, `local_values` — это вход для MPI-сборки.

## 6. Детали реализации

**Файлы:** [all/include/ops_all.hpp](include/ops_all.hpp), [all/src/ops_all.cpp](src/ops_all.cpp).

Поток выполнения `RunImpl`:

```text
1.  rank, size = MPI_Comm_rank, MPI_Comm_size
2.  col_starts[p] = (p * cols_B) / size, p=0..size       // равномерное распределение
3.  j_start, j_end = col_starts[rank], col_starts[rank+1]
4.  ComputeLocalSlab(lhs, rhs, j_start, j_end, ...)      // OMP-параллельно
5.  MPI_Allgather(&local_nnz)                            // сколько ненулей у кого
6.  nnz_displs = prefix-sum(nnz_per_proc)
7.  MPI_Allgatherv(row_indices, values)                  // объединение
8.  local_inner[c] = local_col_ptrs[c] + nnz_displs[rank]
9.  MPI_Allgatherv(local_inner -> out.col_ptrs)
10. out.col_ptrs[cols_B] = total_nnz                     // финальный элемент
```

**Потенциальные узкие места:**

- При `size > cols_B` некоторые ранги получают пустую полосу (`local_cols == 0`).
  В коде есть ветвь `if (local_cols > 0)`, а `local_col_ptrs` инициализируется как `{0}` —
  это корректно работает в `Allgather` (`local_nnz == 0`).
- `MPI_Allgatherv` для `total_nnz` элементов — это `total_nnz * world_size` байт по сети.
  На больших `nnz_C` коммуникация может доминировать.
- `BuildOutputFromColumns`-эквивалент здесь сериальный после `Allgatherv`,
  но он амортизируется по двум коллективам.

## 7. Проверка корректности

Функциональные тесты для `ALL` под обычным запуском (без `mpirun`) **пропускаются**
инфраструктурой курса:

```text
kALL and kMPI tasks are not under mpirun
[ SKIPPED ] ... _all_enabled_TwoByTwoBasic
```

Это нормальное поведение, так как `IsUnderMpirun()` возвращает `false` без `mpirun -n N`.
Проверить функциональную корректность можно запустив тесты под `mpirun -n 1` / `-n 2`.

Производительный тест с `mpirun -n 1/2/4` и `PPC_NUM_THREADS=1/2` проходит — результат на
2000×2000 матрицах совпадает с локально посчитанным эталоном (`ReferenceSpMM`).
Это подтверждает что и MPI-сборка, и сдвиг `col_ptrs` работают корректно.

Согласованность результата между рангами обеспечивается тем, что **все четыре** `Allgatherv`
детерминированно собирают одинаковые данные на каждый ранг.

## 8. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | AMD Ryzen 5 5600 (6 ядер / 12 потоков SMT) |
| RAM | 15 GiB |
| OS | Ubuntu 24.04.2 LTS в Docker (WSL2 ядро 6.6.114.1, хост Windows 11) |
| Compiler | g++ Ubuntu 13.x, `-O3`, `-fopenmp` |
| MPI | OpenMPI 4.1.6 |
| OpenMP | 4.5 (`_OPENMP=201511`) |
| Размер задачи | A, B — банд 2000×2000, ширина ±50 |

Конфигурации запуска:

```bash
# Чистый ALL без mpirun (1 процесс, разное число потоков)
for N in 1 2 4 8; do
  PPC_NUM_THREADS=$N OMP_NUM_THREADS=$N \
    ./build_dilshodov/bin/ppc_perf_tests \
    --gtest_filter="*Dilshodov*all*"
done

# Гибридные конфигурации ranks × threads
PPC_NUM_THREADS=1 OMP_NUM_THREADS=1 mpirun -n 1 \
  ./build_dilshodov/bin/ppc_perf_tests --gtest_filter="*Dilshodov*all*"
PPC_NUM_THREADS=2 OMP_NUM_THREADS=2 mpirun -n 2 \
  ./build_dilshodov/bin/ppc_perf_tests --gtest_filter="*Dilshodov*all*"
PPC_NUM_THREADS=2 OMP_NUM_THREADS=2 mpirun -n 4 \
  ./build_dilshodov/bin/ppc_perf_tests --gtest_filter="*Dilshodov*all*"
```

## 9. Результаты

### 9.1 Одиночный процесс, разное число OMP-потоков (без mpirun)

| Threads | task_run, с | pipeline, с | Speedup vs SEQ@1 | Speedup vs ALL@1×1 | Efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.0258 | 0.0208 | 7.57 | 1.00 | 100% |
| 2 | 0.0160 | 0.0178 | 12.21 | 1.61 | 81% |
| 4 | 0.0082 | 0.0103 | 23.83 | 3.15 | 79% |
| 8 | 0.0065 | 0.0076 | 30.06 | 3.97 | 50% |

### 9.2 Гибридные конфигурации под `mpirun`

| ranks × threads | Total workers | task_run, с | pipeline, с | Speedup vs SEQ@1 | Speedup vs ALL@1×1 | Efficiency |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 × 1 | 1 | 0.0239 | 0.0249 | 8.18 | 1.00 | 100% |
| 2 × 2 | 4 | 0.0092 | 0.0111 | 21.24 | 2.60 | 65% |
| 4 × 2 | 8 | 0.0130 | 0.0153 | 15.03 | 1.84 | 23% |

**Комментарий.**

- Лучшая гибридная конфигурация — `2×2` (2 процесса, 2 потока в каждом): 0.0092 с,
  ускорение ×21 относительно SEQ.
- Конфигурация `4×2` показывает **деградацию** относительно `2×2` (0.0130 vs 0.0092 с)
  при том же общем числе работников (8). Это типичное проявление коммуникационных накладных:
  четыре раунда `Allgatherv` между 4 процессами стоят значимо больше, чем между 2.
- Чистый ALL@8 (8 OMP-потоков, 1 процесс) — 0.0065 с — обгоняет любую гибридную
  конфигурацию. На одной NUMA-ноде (а Ryzen 5 5600 — это одна NUMA-нода) гибрид не даёт
  выигрыша, потому что нет «дальней» памяти, которую MPI разделял бы между процессами.

Цена коммуникации для текущего размера задачи:

- `nnz_C ≈ 400K` ненулей → ~3.2 МБ на `values` (`double`), ~1.6 МБ на `row_indices` (`int`);
- `Allgatherv` для 4 рангов: ≈ 4× этого объёма (каждый получает копию всех) ≈ 19 МБ трафика.
- Время одного `Allgatherv` на shared-memory MPI — порядка миллисекунд,
  что соизмеримо с вычислительной фазой.

## 10. Выводы

ALL-версия достигает максимального ускорения когда вычислительная фаза доминирует над
коммуникационной — это конфигурация с малым числом MPI-процессов и большим числом OMP-потоков
внутри. На текущей одноузловой машине гибридный подход избыточен: чистый OMP/TBB/STL даёт
сопоставимый или лучший результат при меньшей сложности.

Польза гибридной схемы проявилась бы на:

- многоузловом кластере (где MPI — единственный способ задействовать другие узлы);
- задачах с большим объёмом независимых вычислений и малым количеством коммуникаций;
- NUMA-системах, где привязка процесса к NUMA-ноде уменьшает время доступа к памяти.

Ограничения текущей реализации:

1. Каждый ранг хранит полный `lhs` и `rhs` — не масштабируется по памяти.
   Production-вариант делал бы `Bcast(A) + Scatter(B)`.
2. `Allgatherv` отправляет полный результат всем — лишний трафик,
   если потребитель — только один процесс.
3. Сериальная фаза формирования `local_*` и `out.col_ptrs[cols_B] = total_nnz`
   ограничивает закон Амдала, но её доля невелика.
