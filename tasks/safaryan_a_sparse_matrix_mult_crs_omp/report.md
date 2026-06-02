# Sparse Matrix Multiplication

- Student: Safaryan Grigor, group 3823B1PR5
- Technology: OMP
- Variant: 5

## 1. Introduction

This task implements sparse matrix multiplication for double values using the
same CCS-compatible storage model as the existing sequential implementation.

## 2. Parallelization Scheme

The OpenMP version splits result columns between worker threads. Each thread
computes an independent range of columns into a local partial sparse matrix.
After the parallel region, partial matrices are merged in column order to
preserve deterministic output layout.
