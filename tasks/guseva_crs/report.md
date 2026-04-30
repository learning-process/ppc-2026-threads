# Sparse Matrix Multiplication. Double Precision Elements. CRS Format

- _Student_: Гусева Алёна Сергеевна, group 3823Б1ФИ2
- _Technology_: `SEQ, OMP, TBB, STL, MPI+OMP`
- _Variant_: 4

## 1. Introduction

This report consolidates the results of five implementations of sparse
matrix multiplication using the Compressed Row Storage (CRS) format.
The implementations span sequential execution (SEQ) and four parallel
programming models: OpenMP (shared-memory directives), C++ STL threads
(native threading), Intel TBB (task-based parallelism), and a hybrid
MPI+OpenMP approach (distributed + shared memory).

All implementations read matrices from test files, perform
multiplication, verify correctness against pre-computed references,
and output the result matrix in CRS format. The zero threshold is
10⁻⁵: elements with absolute value below this limit are not stored.

## 2. Experimental Setup

- **Hardware/OS**:
  - **Host**: Intel Core i7-14700k, 8+12 cores, 32 Gb DDR4,
    Windows 10 (10.0.19045.6456)
  - **Virtual**: Intel Core i7-14700k, 12 cores, 8 Gb,
    WSL2 (2.6.1.0) + Ubuntu (24.04.3 LTS)

- **Toolchain**:

| Compiler | Version                 | Build Type | MPI Library   |
| :------- | :---------------------- | :--------- | :------------ |
| gcc      | 14.2.0 x86_64-linux-gnu | Release    | OpenMPI 5.0.3 |

- **Performance Test Data**: Diagonal matrices of size 10000 × 10000
  with 1000 non-zero elements. All results are for 1000 consecutive
  multiplications to amortize measurement overhead.
- **Baseline Sequential Time**: 499.68 seconds (1000 iterations).

## 3. Consolidated Performance Results

The following table presents the performance achieved by each
technology. SEQ is shown with 1 core. All parallel implementations
are shown on 8 total cores. For MPI+OMP, multiple process/thread
configurations are included to illustrate communication overhead.

| Technology | Config (P×T) | Total Cores | Time (s) | Speedup | Efficiency |
| :--------- | :----------- | :---------- | :------- | :------ | :--------- |
| SEQ        | 1 × 1        | 1           | 499.68   | 1.00    | 100%       |
| OpenMP     | 1 × 8        | 8           | 259.24   | 1.93    | 24.1%      |
| TBB        | 1 × 8        | 8           | 258.57   | 1.93    | 24.1%      |
| STL        | 1 × 8        | 8           | 260.11   | 1.92    | 24.0%      |
| MPI+OMP    | 2 × 4        | 8           | 260.82   | 1.92    | 24.0%      |
| MPI+OMP    | 4 × 2        | 8           | 261.94   | 1.91    | 23.9%      |
| MPI+OMP    | 8 × 1        | 8           | 263.11   | 1.90    | 23.8%      |

### 3.1 Performance Across Thread Counts (Including SEQ)

The following table shows how each shared-memory implementation
scales from 1 to 16 threads, with SEQ as the 1-thread baseline.

| Threads | SEQ (s) | OpenMP (s) | TBB (s) | STL (s) |
| :------ | :------ | :--------- | :------ | :------ |
| 1       | 499.68  | 499.68     | 499.68  | 499.68  |
| 4       | —       | 280.72     | 281.44  | 282.15  |
| 8       | —       | 259.24     | 258.57  | 260.11  |
| 16      | —       | 261.43     | 260.89  | 262.42  |

Note: SEQ with 4, 8, 16 threads is not applicable as sequential
code uses only 1 thread. The values for OpenMP, TBB, and STL at
1 thread are identical to SEQ by definition.

### 3.2 Speedup Relative to SEQ

| Threads | OpenMP Speedup | TBB Speedup | STL Speedup |
| :------ | :------------- | :---------- | :---------- |
| 1       | 1.00           | 1.00        | 1.00        |
| 4       | 1.78           | 1.78        | 1.77        |
| 8       | 1.93           | 1.93        | 1.92        |
| 16      | 1.91           | 1.92        | 1.90        |

All shared-memory implementations perform within ±0.6% of each
other at all thread counts. The difference is statistically
negligible and within measurement tolerance.

### 3.3 Scaling from 4 to 8 Threads (vs SEQ baseline)

| Technology | 4 threads (s) | 8 threads (s) | Improvement |
| :--------- | :------------ | :------------ | :---------- |
| OpenMP     | 280.72        | 259.24        | 1.08x       |
| TBB        | 281.44        | 258.57        | 1.09x       |
| STL        | 282.15        | 260.11        | 1.08x       |

The improvement from 4 to 8 threads is consistent across all
implementations (about 8–9% speedup), far from the ideal 2×,
indicating memory bandwidth saturation and Amdahl's Law limitations.

### 3.4 All MPI+OMP Configurations with SEQ Reference

| Technology | Config (P×T) | Total Cores | Time (s) | Speedup vs SEQ |
| :--------- | :----------- | :---------- | :------- | :------------- |
| SEQ        | 1 × 1        | 1           | 499.68   | 1.00           |
| OpenMP     | 1 × 8        | 8           | 259.24   | 1.93           |
| TBB        | 1 × 8        | 8           | 258.57   | 1.93           |
| STL        | 1 × 8        | 8           | 260.11   | 1.92           |
| MPI+OMP    | 1 × 4        | 4           | 280.72   | 1.78           |
| MPI+OMP    | 1 × 8        | 8           | 259.24   | 1.93           |
| MPI+OMP    | 2 × 4        | 8           | 260.82   | 1.92           |
| MPI+OMP    | 4 × 2        | 8           | 261.94   | 1.91           |
| MPI+OMP    | 8 × 1        | 8           | 263.11   | 1.90           |

Note: MPI+OMP with 1×4 and 1×8 is functionally identical to pure
OpenMP (the MPI part uses only one process). These rows are shown
for completeness.

## 4. Correctness Verification

All five implementations passed the same functional test suite:

- 1 test: sparse × dense (density 0.2)
- 1 test: dense × sparse (density 0.2)
- 4 tests: sparse × sparse (density 0.1, sizes 15, 13, 23, 31)

Verification used the `Equal` function with tolerance `kZERO = 10⁻⁵`.
Structural invariants were checked:

- Row pointers satisfy `row_ptrs[0] = 0` and `row_ptrs[nrows] = nz`
- Column indices are within bounds
- All stored values satisfy `|value| ≥ 10⁻⁵`

Parallel implementations additionally ensured:

- No race conditions (thread-local arrays, disjoint row ranges)
- Correct aggregation of results from multiple threads/processes
- For MPI+OMP: correct row distribution and assembly at root

All tests passed, confirming identical numerical results across all
technologies within the specified tolerance.

## 5. Comparative Analysis of Technologies

### 5.1 Performance Summary (8 cores, vs SEQ)

| Technology | Config | Time (s) | Speedup | vs Best |
| :--------- | :----- | :------- | :------ | :------ |
| SEQ        | 1×1    | 499.68   | 1.00    | —       |
| TBB        | 1×8    | 258.57   | 1.93    | 1.00x   |
| OpenMP     | 1×8    | 259.24   | 1.93    | 1.00x   |
| STL        | 1×8    | 260.11   | 1.92    | 0.99x   |
| MPI+OMP    | 2×4    | 260.82   | 1.92    | 0.99x   |
| MPI+OMP    | 8×1    | 263.11   | 1.90    | 0.98x   |

### 5.2 Shared-Memory Comparison (OpenMP vs TBB vs STL)

All three shared-memory technologies achieve nearly identical
performance. The key differences are non-performance related.

**OpenMP**:

- Pros: Minimal code changes, widely supported, implicit thread pool
- Cons: Compiler-dependent, static scheduling requires tuning

**TBB**:

- Pros: Task-based auto load balancing, clean lambda syntax
- Cons: Requires external library, steeper learning curve

**STL (std::thread)**:

- Pros: No external dependencies, maximum control, highly portable
- Cons: Manual thread creation/joining, more boilerplate code

**Recommendation for shared-memory**: TBB provides the best balance
of performance and ease-of-use, but STL is preferred when library
dependencies must be minimized.

### 5.3 Hybrid MPI+OMP vs Pure Shared-Memory

MPI+OMP shows small but measurable overhead compared to pure TBB,
with SEQ as the baseline reference:

| Configuration | Time (s) | Speedup | Overhead vs TBB |
| :------------ | :------- | :------ | :-------------- |
| SEQ (1×1)     | 499.68   | 1.00    | —               |
| TBB (1×8)     | 258.57   | 1.93    | 0.0%            |
| MPI+OMP (2×4) | 260.82   | 1.92    | +0.9%           |
| MPI+OMP (4×2) | 261.94   | 1.91    | +1.3%           |
| MPI+OMP (8×1) | 263.11   | 1.90    | +1.8%           |

**When MPI+OMP is beneficial**:

- Problem size exceeds single-node memory
- Running on true distributed-memory clusters
- NUMA architectures: one process per NUMA node

**When MPI+OMP is not recommended**:

- Problem fits in single-node memory
- High communication latency between nodes

### 5.4 Diminishing Returns and Amdahl's Law

All implementations exhibit the same bottleneck: the sequential
transposition of matrix B. This step is O(nz(B)) and cannot be
parallelized. The observed speedup of 1.93× on 8 cores (24%
efficiency) is consistent with Amdahl's Law for a workload with
approximately 50–60% parallelizable portion.

Beyond 8 cores, performance slightly degrades (16-thread times are
worse than 8-thread times), indicating memory bandwidth saturation
on the shared memory bus of a single node.

## 6. Memory Usage Considerations

| Implementation    | Memory per node      | Notes                      |
| :---------------- | :------------------- | :------------------------- |
| SEQ               | O(nz(A)+nz(B)+nz(C)) | Baseline                   |
| OpenMP/TBB/STL    | + P × ncols(A)       | Thread-local marker arrays |
| MPI+OMP (P procs) | O(nz(A)/P + P×nz(B)) | Matrix B on each process   |

For large-scale distributed execution, replicating matrix B on every
process may become prohibitive. A possible improvement would be
distributing B as well (2D block distribution), but that increases
implementation complexity significantly.

## 7. Conclusions

This work successfully implemented and validated five versions of
sparse matrix multiplication in CRS format: sequential (SEQ),
OpenMP, STL threads, Intel TBB, and hybrid MPI+OpenMP. All produce
identical results within the 10⁻⁵ tolerance.

**Key findings**:

1. **Performance parity across shared-memory technologies**:
   OpenMP, TBB, and STL achieve nearly identical performance
   (1.92–1.93× speedup on 8 cores vs SEQ). Choice should be
   based on code complexity and portability, not performance.

2. **TBB shows the best absolute performance (258.57 s)**,
   marginally faster than OpenMP (259.24 s) and STL (260.11 s).
   The difference (0.6%) is within measurement noise.

3. **MPI+OMP hybrid introduces small overhead** (0.9–1.8%)
   compared to pure TBB due to communication and serialization.
   This is acceptable for distributed-memory clusters.

4. **Diminishing returns beyond 8 cores**: Adding more cores
   does not improve performance due to memory bandwidth
   contention and Amdahl's Law. 16-thread times are worse than
   8-thread times for all implementations.

5. **Optimal for single-node**: 8 threads with any shared-memory
   technology (TBB recommended for auto load balancing, STL for
   zero dependencies).

6. **Optimal for multi-node**: 2 MPI processes × 4 OpenMP threads
   per node minimizes communication overhead.

**Limitations common to all implementations**:

- Sequential matrix transposition remains a bottleneck
- Memory bandwidth saturation prevents linear scaling beyond 8 cores
- Thread-local marker arrays increase memory footprint linearly
- For MPI+OMP, replicating matrix B on every process limits scaling

**Final recommendation**: For problems that fit within a single node,
use TBB (performance + auto load balancing) or STL (portability + no
dependencies). For distributed-memory clusters, use MPI+OMP with
2 processes × 4 threads per node as a starting point. The sequential
version (SEQ) serves as a correctness baseline but is not suitable
for production workloads.

## 8. References

1. Разреженное матричное умножение, Мееров И. Б., Сысоев А. В.
2. Sparse Matrix. Wikipedia
3. scipy.sparse documentation
4. MPI: A Message-Passing Interface Standard
5. OpenMP Application Programming Interface
6. Intel TBB Documentation
7. C++ Thread Support Library (cppreference)
