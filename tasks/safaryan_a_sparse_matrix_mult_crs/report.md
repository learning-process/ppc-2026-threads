# Sparse Matrix Multiplication in CCS Format

- Student: Safaryan Grigor
- Group: 3823B1PR5
- Technologies: SEQ, OMP, TBB, STL, ALL
- Variant: 5
- Element type: `double`
- Storage format: CCS, column-compressed storage

## Summary

This task implements sparse matrix multiplication for matrices stored in CCS
format. The matrix data model uses `values`, `row_indices`, and `col_ptrs`.
The task folder name is kept unchanged, although the implemented storage format
is CCS.

## Implementation Paths

- `seq/` contains the sequential baseline.
- `omp/` contains the OpenMP version.
- `tbb/` contains the oneTBB version.
- `stl/` contains the `std::thread` version.
- `all/` contains the combined ALL version.
- `tests/functional/main.cpp` contains shared functional tests.
- `tests/performance/main.cpp` contains shared performance tests.

## Testing

Functional tests run the same fixed CCS cases for all five implementations.
Performance tests use generated sparse matrices and compare the result with a
dense reference multiplication.

## Local Commands

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
$env:PPC_NUM_PROC="2"
$env:PPC_NUM_THREADS="4"
mpiexec -n 2 .\build\bin\ppc_func_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
mpiexec -n 2 .\build\bin\ppc_perf_tests.exe --gtest_filter="*safaryan_a_sparse_matrix_mult_crs*"
```

## Result

All implementations are placed in one task directory and share one CCS data
model, one settings file, and one set of tests.
