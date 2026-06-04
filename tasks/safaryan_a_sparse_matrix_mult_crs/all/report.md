# ALL Implementation

## Purpose

The ALL implementation is the combined version of sparse matrix multiplication.
It works with `double` values and CCS storage using `values`, `row_indices`,
and `col_ptrs`.

## CCS Storage

CCS stores matrix data by columns. The range for column `j` is defined by
`col_ptrs[j] .. col_ptrs[j + 1]`. Values are stored in `values`, and row
numbers are stored in `row_indices`.

## Algorithm

The base multiplication is column-oriented. For each column of `B`, the
implementation accumulates contributions from matching columns of `A` and then
writes the nonzero result entries to CCS arrays.

## Hybrid Scheme

The implementation follows the project ALL pattern. MPI is used for process
synchronization around the measured computation. Each MPI rank builds the full
result locally, which makes validation identical on every process.

Inside each rank, the work is split by result column ranges. OpenMP prepares
independent chunk boundaries, and oneTBB processes chunks in parallel. Each
worker writes to a private partial CCS matrix.

`MPI_Barrier` is used before and after the computation. The first barrier makes
ranks enter the measured part together. The second barrier ensures all ranks
finish before leaving the task stage.

Race conditions are avoided because final CCS arrays are not written during the
parallel region. Partial matrices are merged sequentially in chunk order, so
the result is deterministic.

The number of MPI processes is controlled by `PPC_NUM_PROC`. The number of
threads inside each process is controlled by `PPC_NUM_THREADS`. Local tests were
run with `mpiexec -n 2`.

## Correctness And Performance

Functional tests compare CCS arrays with expected results. Performance tests
run both `pipeline` and `task_run` modes. No benchmark table is included because
timings depend on the local environment.

## Local Commands

```powershell
cmake --build build --config Release --parallel
$env:PPC_NUM_PROC="2"
$env:PPC_NUM_THREADS="4"
mpiexec -n 2 .\build\bin\ppc_func_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
mpiexec -n 2 .\build\bin\ppc_perf_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
```

## Conclusion

The ALL version combines MPI synchronization with thread-level work splitting
and keeps deterministic CCS output.
