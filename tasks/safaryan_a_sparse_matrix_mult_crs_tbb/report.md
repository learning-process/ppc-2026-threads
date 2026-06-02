# Sparse Matrix Multiplication

- Student: Safaryan Grigor, group 3823B1PR5
- Technology: TBB
- Variant: 5

## 1. Introduction

This task implements sparse matrix multiplication for double values using the
same CCS-compatible storage model as the existing sequential implementation.

## 2. Parallelization Scheme

The TBB version splits result columns into independent chunks. Each chunk is
computed in parallel with `tbb::parallel_for` into a local sparse matrix.
After the parallel stage, partial matrices are merged in column order to
preserve deterministic CCS output layout.
