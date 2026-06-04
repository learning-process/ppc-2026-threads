# oneTBB Implementation

## Purpose

The oneTBB implementation multiplies sparse matrices with `double` values in
CCS format. It uses the shared task data model with `values`, `row_indices`,
and `col_ptrs`.

## CCS Storage

In CCS, each column is represented by a range in `col_ptrs`. The values in that
range are stored in `values`, and their row numbers are stored in
`row_indices`. This matches the column-oriented multiplication algorithm.

## Algorithm

Each result column is computed with a temporary accumulator. Nonzero elements of
the current column of `B` select columns of `A`, and their products are added to
the accumulator. The nonzero accumulator entries are stored in CCS form.

## Parallelization

The result columns are split into continuous chunks. The implementation uses
`tbb::parallel_for` over chunk identifiers. Each worker processes only its own
column range and writes to its own partial CCS matrix.

Input matrices are shared read-only data. The final `values`, `row_indices`,
and `col_ptrs` arrays are not modified during the parallel region. After
`tbb::parallel_for` finishes, partial matrices are merged in chunk order. This
prevents data races and preserves deterministic output ordering.

## Correctness And Performance

Functional tests compare the result with expected CCS matrices. Performance
tests run both course modes: `pipeline` for the full task lifecycle and
`task_run` for the main computation stage.

No timing table is included because timings depend on the local environment.

## Local Commands

```powershell
cmake --build build --config Release --parallel
$env:PPC_NUM_THREADS="4"
.\build\bin\ppc_func_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
.\build\bin\ppc_perf_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
```

## Conclusion

The oneTBB version parallelizes independent column chunks, avoids shared writes
inside the parallel region, and produces deterministic CCS output.
