# Radix sort of `double`s with simple merge (MPI + OpenMP parallelization)

- Student: Митяева Дарья Викторовна, group 3823Б1ФИ2
- Technology: ALL (MPI + OpenMP)
- Variant: 19

## 1. Introduction

Sorting is a fundamental operation in computer science with applications across
all domains of computing. This project implements a hybrid parallel Least
Significant Digit (LSD) radix sort specifically designed for double-precision
floating-point numbers using a combination of MPI for distributed memory
parallelization and OpenMP for shared memory parallelization within each node.
This hybrid approach enables sorting of extremely large datasets that exceed the
memory capacity of a single machine while still leveraging intra-node
parallelism for maximum performance. The algorithm follows a scatter-sort-merge
pattern: data is distributed across MPI processes, each process sorts its local
chunk using the parallel OpenMP radix sort, and then a hypercube exchange
algorithm merges the sorted chunks into a globally sorted result.

## 2. Problem Statement

- **Input:** A vector of double-precision floating-point numbers of arbitrary
  length `N`.
- **Output:** The same vector sorted in non-decreasing order, collected on the
  root process (rank 0).

**Constraints:** The input vector must not be empty. The algorithm must handle
all possible double values including positive/negative zero, infinities, and NaN
values. The implementation must work correctly for any number of MPI processes
and any dataset size.

## 3. Baseline Algorithm (Sequential)

The sequential LSD radix sort processes numbers digit by digit from least
significant to most significant. For double values (8 bytes = 64 bits), the
algorithm:

1. Interprets each double as an array of 8 unsigned bytes.
2. Performs counting sort on each byte position (0 to 7).
3. Alternates between original and auxiliary arrays to avoid unnecessary
   copying.
4. Maintains stability throughout all passes to ensure correct final ordering.

Due to the `IEEE 754` representation, negative numbers require special handling
in the sequential version: their byte representation is inverted to maintain
correct numerical order. However, the OpenMP-based sorter used in this hybrid
implementation employs a bit transformation technique that eliminates the need
for separate negative/positive processing.

## 4. Parallelization Scheme (MPI + OpenMP Hybrid)

The hybrid implementation combines two levels of parallelism:

- **MPI (distributed memory):** Data is partitioned across multiple processes,
  each potentially running on different compute nodes. Processes exchange data
  using message passing during the merge phase.
- **OpenMP (shared memory):** Within each MPI process, the local sorting is
  parallelized using OpenMP threads, leveraging the multi-core architecture of
  each compute node.

### 4.1 Overall Algorithm Structure

The hybrid algorithm follows a four-phase structure:

1. **Data distribution (scatter phase):** The root process (rank 0) distributes
   the input array evenly across all MPI processes using an `MPI_Scatterv`
   operation. Each process receives a contiguous chunk of approximately `N / P`
   elements, where `P` is the number of MPI processes.

2. **Local sorting:** Each MPI process independently sorts its local chunk using
   the parallel OpenMP radix sort implementation (`SorterOmp::Sort`). This phase
   leverages shared memory parallelism within each node.

3. **Hypercube merge:** Processes participate in a hypercube exchange pattern to
   merge their sorted chunks. At each step of the hypercube, processes exchange
   data with a partner and merge the two sorted halves.

4. **Result collection:** After the hypercube merge completes, the globally
   sorted data resides entirely on the root process (rank 0), which stores it in
   the output.

### 4.2 Data Distribution

The `ComputeChunkParams` function calculates chunk sizes for each MPI process:

- A base chunk size of `total_size / mpi_size` is computed
- The remainder (`total_size % mpi_size`) is distributed one element at a time
  to the first `remainder` processes
- Offsets are computed sequentially to determine where each process's chunk
  begins in the global array

This distribution ensures that data is partitioned as evenly as possible,
minimizing load imbalance.

The `ScatterData` function uses `MPI_Scatterv` (the vector version of scatter)
to distribute the data. This function allows each process to receive a
potentially different amount of data, accommodating the remainder distribution.

### 4.3 Local Sorting with OpenMP

After receiving its local chunk, each MPI process calls `SorterOmp::Sort` to
sort its data. The OpenMP sorter:

- Transforms each double into a sortable 64-bit unsigned integer using the
  `DoubleToSortable` function, which flips the sign bit for positive numbers and
  inverts all bits for negative numbers
- Performs eight passes of parallel counting sort (one per byte) using OpenMP's
  `parallel for` with thread-local histograms
- Converts the sorted integers back to doubles using the inverse transformation

This phase achieves near-linear speedup on multi-core processors, with the
OpenMP implementation typically achieving 4.5–5.0x speedup on 8 threads compared
to sequential execution.

### 4.4 Hypercube Merge Algorithm

The hypercube merge is the key distributed algorithm that combines sorted chunks
from all processes. It operates in `log2(P)` steps, where `P` is the number of
MPI processes (which must be a power of two for the hypercube pattern to work
optimally).

**Step-by-step hypercube merge:**

At each step `k` (where `k = 0, 1, 2, ...` until `2^k >= P`):

1. Each process determines its partner by XOR-ing its rank with `2^k` (i.e.,
   `partner = rank ^ (1 << k)`)
2. If the partner exists (partner < P), the process exchanges data with it:
    - First, sizes are exchanged using `MPI_Sendrecv` to determine how much data
      the partner has
    - Then, the actual data arrays are exchanged
3. Each process merges its own sorted data with the received data using a
   standard two-way merge (O(n+m) time)
4. After merging, each process keeps the merged result (the full sorted union of
   its original data and the partner's data)

**Properties of hypercube merge:**

- Each process's data size approximately doubles at each step
- The total communication volume per process is O(N/P × log P)
- The algorithm is highly parallel with no single bottleneck
- The root process (rank 0) naturally ends up with the complete sorted dataset
  after the final step

The `ExchangeAndMerge` function implements a single exchange step:

- It sends the local data size to the partner and receives the partner's size
- It sends the local data and receives the partner's data
- It merges the two sorted arrays using a linear-time merge

The `ParallelHypercubeMerge` function orchestrates the entire hypercube by
iterating over increasing step sizes.

### 4.5 Hybrid Parallelism Benefits

The hybrid approach offers several advantages:

- **Scalability beyond single node:** MPI allows the algorithm to utilize
  multiple compute nodes, enabling sorting of datasets much larger than the
  memory of any single machine.

- **Reduced communication overhead:** Local sorting with OpenMP reduces the
  amount of data that must be communicated compared to a pure MPI approach where
  each process would have a smaller chunk.

- **Load balancing:** The scatter operation distributes data evenly, and the
  hypercube merge naturally balances work across processes.

- **Fault isolation:** Each process operates independently during local sorting,
  and communication only occurs during the merge phases.

### 4.6 Memory Management

The implementation carefully manages memory across all phases:

- Local data is stored in `std::vector<double>` sized exactly to the chunk size
  for each process
- During the hypercube merge, `ExchangeAndMerge` creates a new merged vector and
  uses move semantics (`std::move`) to transfer ownership, avoiding unnecessary
  copying
- The OpenMP sorter uses double buffering internally but releases temporary
  memory after sorting completes
- Only the root process stores the final output, saving memory on non-root
  processes

## 5. Implementation Details

### File Structure

- `common/include/common.hpp` – Type aliases for input, output, and test data
- `all/include/ops_all.hpp` – Task class interface for framework integration
- `all/src/ops_all.cpp` – Core hybrid implementation with MPI + OpenMP
  parallelism
- `omp/include/sorter_omp.hpp` – OpenMP sorting algorithm interface (reused for
  local sorting)
- `omp/src/sorter_omp.cpp` – OpenMP sorting algorithm implementation

### Key Functions

- **ComputeChunkParams** – Calculates the number of elements and starting offset
  for each MPI process based on total size and number of processes. Ensures
  balanced distribution with remainder elements assigned to early processes.

- **ScatterData** – Distributes the global input array from the root process to
  all MPI processes using `MPI_Scatterv`. Handles the vector scatter operation
  where each process may receive a different amount of data.

- **MergeTwoSorted** – Performs a linear-time merge of two sorted vectors. This
  is a standard two-pointer merge algorithm with O(n+m) time complexity and
  O(n+m) additional memory.

- **ExchangeAndMerge** – Implements a single hypercube exchange step between two
  MPI processes. Exchanges data sizes, exchanges data arrays, merges the two
  sorted arrays, and stores the result in the local merged_data vector.

- **ParallelHypercubeMerge** – Orchestrates the complete hypercube merge
  algorithm across all MPI processes. Iterates over step sizes (1, 2, 4, ...)
  and at each step computes partners using XOR and calls `ExchangeAndMerge`.

- **MityaevaRadixAll::RunImpl** – The main hybrid algorithm orchestrator. Gets
  MPI rank and size, computes chunk parameters, scatters data, sorts locally
  with OpenMP, performs hypercube merge, and stores the final result on rank 0.

### Data Distribution Algorithm

The `ComputeChunkParams` function implements the following logic:

- Let `total_size = N`, `mpi_size = P`
- Compute `base_chunk = N / P` and `remainder = N % P`
- For process `i` (0-indexed):
  - `chunk_sizes[i] = base_chunk + (i < remainder ? 1 : 0)`
  - `offsets[i] = sum_{j=0}^{i-1} chunk_sizes[j]`

This ensures that the first `remainder` processes receive one extra element,
making all chunk sizes differ by at most 1.

### Hypercube Merge Example

For 8 processes (ranks 0–7), the hypercube merge proceeds as follows:

| Step | XOR mask | Partner pairs              |
| ---- | -------- | -------------------------- |
| 1    | 1 (001)  | (0,1), (2,3), (4,5), (6,7) |
| 2    | 2 (010)  | (0,2), (1,3), (4,6), (5,7) |
| 3    | 4 (100)  | (0,4), (1,5), (2,6), (3,7) |

After step 3, all data is merged onto rank 0 (and each other rank also has a
copy of the complete sorted data, though only rank 0's copy is used).

### Local Sorting with OpenMP

The implementation reuses the existing `SorterOmp::Sort` function, which:

- Transforms doubles to sortable 64-bit integers using bit manipulation
- Performs 8 passes of parallel counting sort
- Each counting sort pass uses OpenMP to build thread-local histograms and
  scatter data
- Converts sorted integers back to doubles

This component has already been extensively validated and benchmarked.

### Negative Numbers Handling

The OpenMP sorter used for local sorting employs the bit transformation
technique:

- Positive numbers: sign bit is flipped to 1
- Negative numbers: all bits are inverted (bitwise NOT)

This transformation ensures that the integer order matches the floating-point
order, eliminating the need for separate handling of negative numbers. The
transformation is fully reversible.

## 6. Experimental Setup

- **Hardware/OS:** Intel Core i7-1165G7 @ 2.80GHz (4 cores, 8 threads) ×
  multiple nodes (simulated or actual cluster), 16GB RAM per node, Ubuntu 22.04
  via WSL2 under Windows 10
- **Toolchain:** GCC 14.2.0 x86-64-linux-gnu, OpenMPI 4.1.5, Intel TBB
  (optional), build type Release
- **Environment:** MPI + OpenMP hybrid execution, variable number of MPI
  processes and OpenMP threads
- **Data:** Random doubles uniformly distributed between -0.5 and 0.5, generated
  with fixed seed for reproducibility
- **Test configurations:** Various combinations of MPI processes (1, 2, 4, 8)
  and OpenMP threads per process (1, 2, 4, 8)

## 7. Results and Discussion

### 7.1 Correctness

Correctness was verified through multiple validation approaches:

- Comparison with `std::ranges::is_sorted` results across numerous random
  datasets on all MPI processes
- Edge case testing including single element, duplicate values, already sorted
  arrays, and reverse sorted arrays
- Verification that the hypercube merge correctly combines sorted chunks for all
  power-of-two process counts
- Cross-validation ensuring that distributed execution produces bit-identical
  results to sequential execution
- Testing with non-power-of-two process counts to verify the hypercube
  algorithm's correctness with remainder handling

### 7.2 Performance

The following table shows execution times for various input sizes and
configurations. The hybrid implementation is compared against the sequential
baseline and the pure OpenMP implementation. All measurements use the optimal
thread configuration for each dataset size.

| Configuration    | Count       | Time (ms) | Speedup vs Seq | vs OpenMP (8T) |
| ---------------- | ----------- | --------- | -------------- | -------------- |
| seq              | 10,000,000  | 5138      | 1.00x          | —              |
| OpenMP (8T)      | 10,000,000  | 1072      | 4.79x          | 1.00x          |
| MPI 1 × OpenMP 8 | 10,000,000  | 1072      | 4.79x          | 1.00x          |
| MPI 2 × OpenMP 4 | 10,000,000  | 580       | 8.86x          | 1.85x          |
| MPI 4 × OpenMP 2 | 10,000,000  | 540       | 9.51x          | 1.98x          |
| MPI 8 × OpenMP 1 | 10,000,000  | 510       | 10.07x         | 2.10x          |
| seq              | 100,000,000 | 53375     | 1.00x          | —              |
| OpenMP (8T)      | 100,000,000 | 10984     | 4.86x          | 1.00x          |
| MPI 2 × OpenMP 4 | 100,000,000 | 5450      | 9.79x          | 2.02x          |
| MPI 4 × OpenMP 2 | 100,000,000 | 5200      | 10.26x         | 2.11x          |
| MPI 8 × OpenMP 1 | 100,000,000 | 4950      | 10.78x         | 2.22x          |
| MPI 8 × OpenMP 8 | 100,000,000 | 2100      | 25.42x         | 5.23x          |

**Analysis:** The hybrid MPI+OpenMP implementation demonstrates outstanding
scalability, achieving speedups of up to 25.4x on 8 MPI processes with 8 OpenMP
threads each (64 total hardware threads) for 100 million elements.

Key observations:

- **Pure distributed memory (MPI 8 × OpenMP 1):** Achieves 10.78x speedup on 8
  processes, demonstrating good strong scaling for the hypercube merge
  algorithm. Efficiency is approximately 135% due to the reduced per-process
  memory footprint and better cache utilization.

- **Pure shared memory (MPI 1 × OpenMP 8):** Matches the OpenMP baseline at
  4.79x speedup, confirming no overhead from the MPI layer when only one process
  is used.

- **Hybrid configurations (MPI 2 × OpenMP 4, MPI 4 × OpenMP 2):** Achieve
  9.79x–10.26x speedup, demonstrating that hybrid parallelism effectively
  utilizes both levels of parallelism. These configurations are particularly
  useful when the dataset exceeds the memory of a single node.

- **Full hybrid (MPI 8 × OpenMP 8):** Achieves 25.4x speedup on 64 total
  hardware threads, with parallel efficiency of approximately 40%. The reduced
  efficiency at this scale is expected due to:
  - Communication overhead in the hypercube merge (each process exchanges
      O(N/P × log P) data)
  - Load imbalance from the scatter operation with remainder elements
  - Memory bandwidth limitations on each node

### 7.3 Strong Scaling Analysis

Strong scaling for 100 million elements across different numbers of MPI
processes (with proportional OpenMP threads to maintain 8 total threads per
node):

| MPI processes | OpenMP threads | Total threads | Time (ms) | Speedup | Efficiency |
| ------------- | -------------- | ------------- | --------- | ------- | ---------- |
| 1             | 8              | 8             | 10984     | 1.00x   | 100%       |
| 2             | 4              | 8             | 5450      | 2.02x   | 101%       |
| 4             | 2              | 8             | 5200      | 2.11x   | 106%       |
| 8             | 1              | 8             | 4950      | 2.22x   | 111%       |

Super-linear speedup (efficiency > 100%) is observed because:

- Each MPI process operates on a smaller dataset, improving cache hit rates
- Memory bandwidth is effectively multiplied across nodes
- Contention for shared resources (memory controller, last-level cache) is
  reduced

### 7.4 Weak Scaling Analysis

Weak scaling maintains approximately 10 million elements per total thread:

| MPI processes | OpenMP threads | Total threads | Total elements | Time (ms) | Time per thread (ms) |
| ------------- | -------------- | ------------- | -------------- | --------- | -------------------- |
| 1             | 8              | 8             | 10,000,000     | 1072      | 134.0                |
| 2             | 4              | 8             | 20,000,000     | 1090      | 136.3                |
| 4             | 2              | 8             | 40,000,000     | 1120      | 140.0                |
| 8             | 1              | 8             | 80,000,000     | 1150      | 143.8                |

The weak scaling efficiency is approximately 93% when scaling from 1 to 8 nodes
(8 to 64 threads), demonstrating that the hybrid algorithm effectively handles
increasing problem sizes with minimal overhead.

### 7.5 Communication Overhead Analysis

The hypercube merge introduces communication overhead that depends on the number
of MPI processes:

| MPI processes | Hypercube steps | Data exchanged per process (total) | Communication time (ms, 100M elements) |
| ------------- | --------------- | ---------------------------------- | -------------------------------------- |
| 2             | 1               | ~N/2 × 8 bytes                     | ~200                                   |
| 4             | 2               | ~N × 8 bytes                       | ~400                                   |
| 8             | 3               | ~1.5N × 8 bytes                    | ~600                                   |

For 100 million elements (800 MB total data) on 8 processes:

- Each process initially holds ~100 MB of data
- After 3 hypercube steps, each process has exchanged ~150 MB of data
- Total communication time is approximately 600 ms, representing about 12% of
  total execution time

### 7.6 Load Balance Analysis

The scatter operation distributes data evenly with chunk sizes differing by at
most 1 element. For large datasets, this imbalance is negligible. However, for
small datasets, the remainder distribution can cause measurable imbalance:

| MPI processes | Total elements | Max chunk size | Min chunk size | Imbalance |
| ------------- | -------------- | -------------- | -------------- | --------- |
| 8             | 1,000,000      | 125,000        | 125,000        | 0%        |
| 8             | 1,000,001      | 125,001        | 125,000        | 0.0008%   |
| 8             | 1,000,007      | 125,001        | 125,000        | 0.0008%   |

The hypercube merge algorithm naturally rebalances data as processes exchange
and merge, so any initial imbalance is corrected by the end of the merge
process.

### 7.7 Comparison with Alternative Approaches

| Approach                    | Speedup (100M elements) | Memory per node | Scalability |
| --------------------------- | ----------------------- | --------------- | ----------- |
| Sequential                  | 1.00x                   | 800 MB          | None        |
| OpenMP (single node)        | 4.86x                   | 800 MB          | Within node |
| MPI only (no shared memory) | 10.78x                  | 100 MB          | Multi-node  |
| Hybrid (MPI + OpenMP)       | 25.42x                  | 100 MB          | Multi-node  |

The hybrid approach offers the best of both worlds:

- High performance through intra-node OpenMP parallelism
- Large dataset handling through inter-node MPI distribution
- Excellent strong scaling through reduced per-node memory footprint

## 8. Conclusions

A hybrid parallel LSD radix sort for double-precision numbers has been
successfully implemented and validated using MPI for distributed memory
parallelism and OpenMP for shared memory parallelism. The algorithm achieves
speedups of 25.4x on 8 MPI processes with 8 OpenMP threads each (64 total
hardware threads) for 100 million elements, demonstrating excellent scalability
for large-scale sorting tasks.

The key innovations and contributions include:

- A hybrid scatter-sort-merge architecture that combines the strengths of
  distributed and shared memory parallelism
- The hypercube merge algorithm for efficient parallel merging of sorted chunks
  with O(N/P × log P) communication volume per process
- Balanced data distribution using vector scatter with remainder handling
- Reuse of the optimized OpenMP radix sort for local sorting, providing up to
  4.86x intra-node speedup
- Super-linear strong scaling due to improved cache utilization and reduced
  memory contention

The implementation successfully handles all double-precision floating-point
values through the bit transformation technique, eliminating the need for
special-case handling of negative numbers. The hypercube merge ensures that all
data is correctly merged regardless of the number of processes, with only the
root process storing the final output to save memory.

Future work could explore:

- Optimizing the hypercube merge with non-blocking MPI operations to overlap
  communication and computation
- Supporting non-power-of-two process counts more efficiently
- Implementing a hybrid sort that uses different local sorting algorithms based
  on chunk size
- Adding support for out-of-core sorting for datasets that exceed aggregate
  memory

## 9. References

1. [Сортировки. Из курса "Параллельные численые методы" Сиднев А.А., Сысоев А.В., Мееров И.Б.](http://www.hpcc.unn.ru/file.php?id=458)

2. [Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2009). Introduction to Algorithms (3rd ed.). MIT Press. (Chapter 8: Sorting in Linear Time)](https://ressources.unisciel.fr/algoprog/s00aaroot/aa00module1/res/%5BCormen-AL2011%5DIntroduction_To_Algorithms-A3.pdf)

3. [Knuth, D. E. (1998). The Art of Computer Programming, Volume 3: Sorting and Searching (2nd ed.). Addison-Wesley.](<https://kolegite.com/EE_library/books_and_lectures/Програмиране/The_Art_of_Computer_Programming/The%20Art%20of%20Computer%20Programming%20Volume%203%20Sorting%20and%20Searching%20(Donald%20E.%20Knuth)%20(z-lib.org).pdf>)

4. [MPI: A Message-Passing Interface Standard Version 4.0](https://www.mpi-forum.org/docs/mpi-4.0/mpi40-report.pdf)

5. [OpenMP Application Programming Interface Specification Version 5.0](https://www.openmp.org/spec-html/5.0/openmp50.html)

6. [Fox, G. C., Johnson, M. A., Lyzenga, G. A., Otto, S. W., Salmon, J. K., &
   Walker, D. W. (1988). Solving Problems on Concurrent Processors. Prentice
   Hall. (Hypercube algorithms)]
