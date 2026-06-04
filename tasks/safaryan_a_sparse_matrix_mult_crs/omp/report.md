# OpenMP Implementation

## Purpose

The OpenMP implementation multiplies sparse matrices with `double` values in
CCS format. It keeps the common task data model based on `values`,
`row_indices`, and `col_ptrs`.

## CCS Storage

CCS stores nonzero values by columns. For a column `j`, its elements are stored
in the range `col_ptrs[j] .. col_ptrs[j + 1]`. The `values` array stores the
nonzero values, and `row_indices` stores the row for each value.

## Algorithm

The result is built by columns. For each result column, a temporary dense
accumulator of length `A.rows` is used. Contributions are added from matching
columns of `A` and nonzero values of the current column of `B`. Nonzero
accumulator entries are then written back to CCS arrays.

## Parallelization

The column range of the result matrix is split into independent chunks. OpenMP
parallelizes the loop over these chunks. Input matrices are shared and read
only. Each worker owns a private accumulator and a private partial CCS matrix.

Race conditions are avoided because workers never write to the final result
during the parallel region. After the parallel region, partial matrices are
merged sequentially in increasing chunk order. This keeps the CCS output
deterministic.

The number of workers is controlled by the course helper through
`PPC_NUM_THREADS`. The `OMP_NUM_THREADS` variable may also be used locally, but
the tests primarily use `PPC_NUM_THREADS`.

## Correctness And Performance

Functional tests compare matrix dimensions, `values`, `row_indices`, and
`col_ptrs` with expected CCS results. Performance tests run both `pipeline` and
`task_run` modes and compare the result with a dense reference multiplication.

No benchmark table is included because timings depend on the local environment.

## Local Commands

```powershell
cmake --build build --config Release --parallel
$env:PPC_NUM_THREADS="4"
.\build\bin\ppc_func_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
.\build\bin\ppc_perf_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
```

## Conclusion

The OpenMP version uses independent column chunks and local partial buffers to
avoid races while preserving deterministic CCS output.
