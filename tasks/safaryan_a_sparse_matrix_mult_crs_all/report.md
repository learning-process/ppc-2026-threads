# Sparse Matrix Multiplication

- Student: Safaryan Grigor, group 3823B1PR5
- Technology: ALL
- Variant: 5

## 1. Introduction

This task implements sparse matrix multiplication for double values using the
same CCS-compatible storage model as the existing sequential implementation.

## 2. Parallelization Scheme

The ALL version follows the project pattern for `kALL` tasks. MPI barriers
synchronize ranks around the measured computation, OpenMP prepares independent
column chunks, and TBB computes those chunks in parallel. Each rank builds the
full CCS result so functional and performance checks can validate the same
output on every process.
