# Умножение разреженных матриц (CCS, double) — SEQ

- Student: Дилшодов Адхам Умидович, группа 3823Б1ПР4
- Technology: SEQ (последовательная реализация)
- Variant: 5

## 1. Контекст

Данный отчёт описывает последовательную реализацию умножения двух разреженных матриц `C = A · B`,
хранящихся в столбцовом разреженном формате CCS (Compressed Column Storage).
SEQ-версия играет роль эталона корректности и нижней границы по времени:
именно она используется в качестве baseline во всех остальных backend-отчётах при подсчёте speedup.

## 2. Постановка задачи

**Вход:** две разреженные матрицы `A` (`m × k`) и `B` (`k × n`) в формате CCS:

```cpp
struct SparseMatrixCCS {
  int rows_count;                  // число строк
  int cols_count;                  // число столбцов
  int non_zeros;                   // число ненулевых элементов
  std::vector<int>    col_ptrs;    // длина cols_count + 1
  std::vector<int>    row_indices; // длина non_zeros
  std::vector<double> values;      // длина non_zeros
};
```

**Выход:** разреженная матрица `C` (`m × n`) в том же CCS-формате.

**Ограничения:**

- `A.cols_count == B.rows_count` (согласованность размеров);
- `col_ptrs[0] == 0`, `col_ptrs[cols_count] == values.size()` (структурная инвариантность);
- `row_indices` строго возрастают внутри каждого столбца;
- все индексы строк лежат в `[0, rows_count)`.

**Критерий корректности:** при умножении на эталон, посчитанный из плотного представления,
выходная матрица должна совпадать по `col_ptrs` и `row_indices` побитно
и иметь `|C_ij - C_ij^ref| < 10⁻⁸` по значениям.

## 3. Базовый алгоритм

Используется классическая схема умножения «столбец на столбец»: каждый столбец `j` результата `C`
получается как линейная комбинация столбцов `A`, индексированных ненулями `j`-го столбца `B`.

```text
for j in 0..n-1:
    acc = пустой словарь {row -> value}
    for каждое ненулевое (k, B[k][j]) в col_j матрицы B:
        for каждое ненулевое (i, A[i][k]) в col_k матрицы A:
            acc[i] += A[i][k] * B[k][j]
    выгрузить acc как j-й столбец C (сортированный по i, с отбрасыванием численных нулей)
```

В реализации `acc` — это `std::map<int, double>`, что автоматически даёт сортировку по `row_index`.

**Асимптотика:**

- Время: `O(Σ_j Σ_{k∈B_j} nnz_A_k · log m)` за счёт `std::map`.
  В худшем случае — `O(flop · log m)`, где `flop` — количество скалярных произведений.
- Память: `O(m)` на словарь-аккумулятор + `O(nnz_C)` на результат.

**Инварианты, проверяемые на входе:**

- `rows_count > 0 && cols_count > 0`;
- `col_ptrs.size() == cols_count + 1`, `col_ptrs.front() == 0`, `col_ptrs.back() == values.size()`;
- `row_indices` строго возрастают в каждом столбце, индексы внутри `[0, rows_count)`.

## 4. Детали реализации

**Файлы:** [seq/include/ops_seq.hpp](include/ops_seq.hpp), [seq/src/ops_seq.cpp](src/ops_seq.cpp).

**`ValidationImpl`** проверяет три группы свойств обеих входных матриц через функции
`HasValidDimensions`, `HasValidContainers`, `HasValidColumnOrdering`,
после чего сверяет согласованность размеров `A.cols == B.rows`.

**`PreProcessingImpl`** обнуляет выходную матрицу (`GetOutput() = SparseMatrixCCS{}`),
чтобы старые данные не утекли между запусками.

**`RunImpl`** — основной цикл по столбцам B:

```cpp
// File: seq/src/ops_seq.cpp
for (int col_b = 0; col_b < matrix_b.cols_count; ++col_b) {
  std::map<int, double> accumulator;
  for (int idx_b = matrix_b.col_ptrs[col_b]; idx_b < matrix_b.col_ptrs[col_b + 1]; ++idx_b) {
    const int pivot_row = matrix_b.row_indices[idx_b];
    const double pivot_value = matrix_b.values[idx_b];
    for (int idx_a = matrix_a.col_ptrs[pivot_row]; idx_a < matrix_a.col_ptrs[pivot_row + 1]; ++idx_a) {
      const int result_row = matrix_a.row_indices[idx_a];
      accumulator[result_row] += matrix_a.values[idx_a] * pivot_value;
    }
  }
  for (const auto &[row, value] : accumulator) {
    if (std::abs(value) > kEps) {
      matrix_c.row_indices.push_back(row);
      matrix_c.values.push_back(value);
    }
  }
  matrix_c.col_ptrs[col_b + 1] = static_cast<int>(matrix_c.values.size());
}
```

Поскольку `std::map` упорядочен по ключу, обход в `for (const auto &[row, value] : accumulator)`
уже даёт строки в возрастающем порядке — это удовлетворяет инварианту CCS на выходе.
Численный нуль (`|value| < 10⁻¹⁰`) отбрасывается, чтобы результат оставался разреженным.

**`PostProcessingImpl`** просто синхронизирует поле `non_zeros` с фактическим размером `values`.

## 5. Проверка корректности

В функциональных тестах ([tests/functional/main.cpp](../tests/functional/main.cpp)) три параметра:
`TwoByTwoBasic`, `ThreeByThreeSparse`, `RectangularCheck`. Каждый принимает плотную матрицу,
конвертируется в CCS через `DenseToCCS`, а эталон считается обычным плотным умножением `DenseMul`
и тоже переводится в CCS. Сравнение через `CompareCCS` строгое: совпадают размеры,
`col_ptrs`, `row_indices`; значения сверяются с точностью `1e-8`.

Все три теста для SEQ проходят (см. лог: `seq_enabled_TwoByTwoBasic`,
`_ThreeByThreeSparse`, `_RectangularCheck` — `[ OK ] (0 ms)`).

## 6. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | AMD Ryzen 5 5600 (6 ядер / 12 потоков) |
| RAM | 15 GiB |
| OS | Ubuntu 24.04.2 LTS внутри Docker, ядро WSL2 6.6.114.1 (хост — Windows 11) |
| Compiler | g++ из дистрибутива Ubuntu 24.04 (GCC 13.x), `-O3` (CMAKE_BUILD_TYPE=Release) |
| OpenMP | 4.5 (`_OPENMP=201511`) |
| MPI | OpenMPI 4.1.6 |
| CMake | 3.28.3 |
| Размер задачи (perf) | A,B — банд-матрицы 2000×2000, ширина ±50, ~200K ненулей, seed=42 |

**Команда запуска SEQ:**

```bash
cmake -S . -B build_dilshodov \
  -DPPC_IMPLEMENTATIONS="all;mpi;omp;seq;stl;tbb" \
  -DPPC_TASKS=dilshodov_a_spmm_double_css
cmake --build build_dilshodov --config Release --target ppc_perf_tests -j 12
./build_dilshodov/bin/ppc_perf_tests --gtest_color=no \
  --gtest_filter="*Dilshodov*seq*"
```

## 7. Результаты

SEQ не масштабируется по потокам — `PPC_NUM_THREADS` для неё не имеет значения,
время определяется только характером алгоритма и шумом измерений.
Ниже среднее по нескольким прогонам инфраструктуры `BaseRunPerfTests`:

| Конфигурация | task_run, с | pipeline, с |
| --- | ---: | ---: |
| SEQ, p=1 | 0.1954 | 0.1999 |
| SEQ, p=2 (формально, тот же бинарь) | 0.1974 | 0.1876 |
| SEQ, p=4 | 0.2206 | 0.2090 |
| SEQ, p=8 | 0.1952 | 0.1925 |

Главный наблюдаемый эффект — самые дорогие операции это **вставка в `std::map`**
(для каждой пары `(i, k)` происходит лог-поиск в красно-чёрном дереве)
и **повторное создание/разрушение** дерева для каждой колонки B.
Именно из-за `std::map` ускорение последующих параллельных реализаций относительно SEQ
значительно превышает наивную оценку «по числу потоков»: параллельные backend-ы заменили
словарь на массив-аккумулятор с маркерами, что само по себе даёт примерно ×8 ускорения
даже на одном потоке.

## 8. Выводы

SEQ-реализация выбрана за простоту и читаемость, а не за скорость. Она служит:

1. эталоном корректности — все backend-ы сравниваются с её выходом;
2. формальным baseline `T_seq` для подсчёта speedup и эффективности в OMP, TBB, STL и ALL.

При обсуждении speedup в сводном отчёте важно отделять *алгоритмический выигрыш*
(sparse accumulator вместо `std::map`) от *параллельного выигрыша*
(распределение колонок по потокам/процессам). Поэтому в локальных отчётах OMP/TBB/STL/ALL
приводятся две оценки ускорения: относительно SEQ и относительно того же backend-а
при одном потоке.
