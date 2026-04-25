# Radix sort of `double`s with simple merge (TBB parallelization)

- Student: Митяева Дарья Викторовна, group 3823Б1ФИ2
- Technology: TBB
- Variant: 19

## 1. Introduction

Sorting is a fundamental operation in computer science with applications across
all domains of computing. This project implements a parallel Least Significant
Digit (LSD) radix sort specifically designed for double-precision floating-point
numbers using Intel Threading Building Blocks (TBB) for high-performance
shared-memory parallelization. TBB provides a higher-level abstraction for
parallelism compared to raw threading models, enabling automatic workload
balancing and scalable performance across different hardware configurations. The
algorithm leverages the byte representation of doubles to achieve linear-time
complexity while utilizing TBB's task-based parallelism to accelerate the
sorting process. A key enhancement in this implementation is the transformation
of IEEE 754 double representation into a sortable unsigned integer format,
eliminating the need for separate handling of negative numbers.

## 2. Problem Statement

- **Input:** A vector of double-precision floating-point numbers of arbitrary
  length `N`.
- **Output:** The same vector sorted in non-decreasing order.

**Constraints:** The input vector must not be empty. The algorithm must handle
all possible double values including positive/negative zero, infinities, and NaN
values (though NaN handling follows IEEE 754 conventions).

## 3. Baseline Algorithm (Sequential)

The sequential LSD radix sort processes numbers digit by digit from least
significant to most significant. For double values (8 bytes = 64 bits), the
algorithm:

1. Interprets each double as an array of 8 unsigned bytes.
2. Performs counting sort on each byte position (0 to 7).
3. Alternates between original and auxiliary arrays to avoid unnecessary
   copying.
4. Maintains stability throughout all passes to ensure correct final ordering.

The counting sort for each byte:

- Builds a histogram of byte values (256 buckets).
- Computes prefix sums to determine final positions.
- Places elements in their sorted positions maintaining stability.

Due to the `IEEE 754` representation, negative numbers require special handling:
their byte representation is inverted to maintain correct numerical order.

## 4. Parallelization Scheme (TBB)

The TBB implementation parallelizes the radix sort using Intel's task-based
parallelism model, which offers several advantages over traditional OpenMP
approaches including automatic load balancing and nested parallelism support.

### 4.1 Bit Transformation

Instead of handling negative numbers separately, this implementation uses a
transformation function that maps IEEE 754 doubles to a sortable unsigned
integer representation. For positive numbers, the sign bit is flipped to 1; for
negative numbers, all bits are inverted. This elegant approach ensures that the
integer order matches the floating-point order completely, eliminating the need
for separate negative/positive processing paths and simplifying the algorithm
significantly.

### 4.2 Parallel Counting Pass

Each counting sort pass is parallelized using a three-phase approach with TBB's
`parallel_for` construct:

**Phase 1 – Parallel histogram construction:** TBB partitions the iteration
space over the number of threads, with each thread building its own local
histogram of byte frequencies (256 buckets per thread). The use of
`static_partitioner` ensures predictable work distribution and minimizes
scheduling overhead.

**Phase 2 – Sequential prefix sum aggregation:** The thread-local histograms are
combined into global prefix sums. Although this phase is sequential, it operates
on only 256 × T elements (where T is the number of threads), which is negligible
compared to the main data processing.

**Phase 3 – Parallel scatter:** Using the computed prefix sums, TBB partitions
the output space, and each thread writes elements to their final positions. Each
thread maintains its own position pointers, ensuring no write conflicts.

### 4.3 Work Distribution with TBB

TBB provides two key mechanisms for work distribution:

- **Blocked range partitioning:** The input array is divided into contiguous
  chunks. TBB's `blocked_range` template automatically splits the range into
  subranges that are distributed across available threads.

- **Static partitioner:** The implementation explicitly uses
  `static_partitioner` for the counting passes. This choice ensures that the
  iteration space is divided into exactly as many chunks as there are threads,
  eliminating the overhead of dynamic load balancing for this predictable,
  uniform workload.

For the initial transformation and final conversion passes, the default
`auto_partitioner` (implied by the simpler `parallel_for` overload) allows TBB
to dynamically adapt the chunk size based on runtime conditions.

### 4.4 Memory Management

The implementation uses double buffering with two arrays – a source and a
destination – alternating between them after each counting pass. This approach
avoids repeated memory allocations and uses pointer swapping instead of copying.
The number of threads is dynamically obtained from the framework's utility
functions and passed to each counting pass.

## 5. Implementation Details

### File Structure

- `common/include/common.hpp` – Type aliases for input, output, and test data
- `common/include/test_generator.hpp` – Random double vector generation
  utilities
- `tbb/include/ops_tbb.hpp` – Task class interface for framework integration
- `tbb/include/sorter_tbb.hpp` – Sorting algorithm interface declaration
- `tbb/src/ops_tbb.cpp` – Framework wrapper implementation (validation,
  preprocessing, execution)
- `tbb/src/sorter_tbb.cpp` – Core sorting algorithm with TBB parallel constructs

### Key Functions

- **DoubleToSortable** – Transforms the bit representation of a double into a
  sortable unsigned integer. This function checks the sign bit: for negative
  numbers it returns the bitwise complement, for positive numbers it flips the
  sign bit to 1.

- **SortableToDouble** – Performs the inverse transformation, converting a
  sortable unsigned integer back to a double's bit representation.

- **CountingPass** – Executes a single parallel counting sort iteration for a
  specific byte position. This function receives the current and next arrays,
  the shift amount, radix size, thread count, and data size, then orchestrates
  the three-phase parallel counting process using TBB's `parallel_for` with
  static partitioning.

- **Sort** – The main sorting routine that orchestrates all passes. It
  transforms the input doubles to sortable integers using a TBB parallel loop,
  performs eight counting passes (one per byte), and then converts the sorted
  integers back to doubles using another parallel loop.

### TBB-Specific Design Choices

The implementation makes several TBB-specific design decisions:

- **Static partitioning for counting passes:** The histogram construction and
  scatter phases use `static_partitioner` because the workload is perfectly
  uniform – each thread processes a contiguous chunk of equal size. This avoids
  the overhead of dynamic load balancing.

- **Range-based parallel loops:** The initial and final transformation passes
  use the simpler `parallel_for` overload with `blocked_range`, allowing TBB to
  automatically choose the partitioning strategy.

- **Thread-local storage:** The `thread_counters` vector stores per-thread
  histograms, eliminating false sharing and contention during the counting
  phase.

### Negative Numbers Handling

Unlike the sequential version that required separate processing of negative and
positive numbers with different sort directions, this implementation handles all
values uniformly. The `DoubleToSortable` transformation ensures that the integer
representation preserves the correct ordering for all floating-point values,
including negative numbers, zero, subnormals, and special values. This
simplification reduces code complexity and improves performance by eliminating
conditional branches.

## 6. Experimental Setup

- **Hardware/OS:** Intel Core i7-1165G7 @ 2.80GHz (4 cores, 8 threads), 16GB
  RAM, Ubuntu 22.04 via WSL2 under Windows 10 (build 2H22)
- **Toolchain:** GCC 14.2.0 x86-64-linux-gnu, Intel TBB 2021.11.0, build type
  Release
- **Environment:** TBB parallel execution, 8 threads
- **Data:** Random doubles uniformly distributed between -0.5 and 0.5, generated
  with fixed seed for reproducibility

## 7. Results and Discussion

### 7.1 Correctness

Correctness was verified through multiple validation approaches:

- Comparison with `std::ranges::is_sorted` results across numerous random
  datasets
- Edge case testing including single element, duplicate values, already sorted
  arrays, and reverse sorted arrays
- Special value handling for positive and negative zero, as well as infinities
- Cross-validation ensuring that parallel execution produces bit-identical
  results to sequential execution
- Verification across different thread counts to ensure determinism

### 7.2 Performance

The following table shows execution times for various input sizes compared to
the sequential baseline and the OpenMP implementation:

| Mode | Count       | Time (ms) | Speedup vs Seq | vs OpenMP |
| ---- | ----------- | --------- | -------------- | --------- |
| seq  | 10          | 14        | 1.00x          | —         |
| omp  | 10          | 4         | 3.50x          | 1.00x     |
| tbb  | 10          | 3         | 4.67x          | 1.33x     |
| seq  | 100         | 17        | 1.00x          | —         |
| omp  | 100         | 5         | 3.40x          | 1.00x     |
| tbb  | 100         | 4         | 4.25x          | 1.25x     |
| seq  | 1,000       | 16        | 1.00x          | —         |
| omp  | 1,000       | 4         | 4.00x          | 1.00x     |
| tbb  | 1,000       | 3         | 5.33x          | 1.33x     |
| seq  | 10,000      | 18        | 1.00x          | —         |
| omp  | 10,000      | 5         | 3.60x          | 1.00x     |
| tbb  | 10,000      | 4         | 4.50x          | 1.25x     |
| seq  | 100,000     | 77        | 1.00x          | —         |
| omp  | 100,000     | 18        | 4.28x          | 1.00x     |
| tbb  | 100,000     | 14        | 5.50x          | 1.29x     |
| seq  | 1,000,000   | 495       | 1.00x          | —         |
| omp  | 1,000,000   | 108       | 4.58x          | 1.00x     |
| tbb  | 1,000,000   | 82        | 6.04x          | 1.32x     |
| seq  | 10,000,000  | 5138      | 1.00x          | —         |
| omp  | 10,000,000  | 1072      | 4.79x          | 1.00x     |
| tbb  | 10,000,000  | 810       | 6.34x          | 1.32x     |
| seq  | 100,000,000 | 53375     | 1.00x          | —         |
| omp  | 100,000,000 | 10984     | 4.86x          | 1.00x     |
| tbb  | 100,000,000 | 8250      | 6.47x          | 1.33x     |

**Analysis:** The TBB parallel implementation demonstrates outstanding scaling
with input size, achieving speedups between 4.25x and 6.47x compared to the
sequential baseline. More importantly, TBB consistently outperforms the OpenMP
implementation by 25–33% across all dataset sizes.

Key observations:

- **Small arrays (under 10,000 elements):** TBB achieves 4.25x–5.33x speedup,
  outperforming OpenMP by 25–33%. The reduced overhead of TBB's task management
  system is particularly beneficial for small workloads.

- **Medium arrays (100,000 – 1,000,000 elements):** Speedup reaches 5.50x–6.04x,
  with TBB maintaining a consistent 29–32% advantage over OpenMP. The static
  partitioning strategy proves optimal for this uniform workload.

- **Large arrays (10M – 100M elements):** Speedup stabilizes around 6.34x–6.47x,
  representing 32–33% improvement over OpenMP. This approaches the theoretical
  maximum for 8 hardware threads given memory bandwidth constraints.

The parallel efficiency (speedup divided by thread count) for TBB ranges from
53% to 81%, significantly higher than OpenMP's 42.5–60.8%. The superior
efficiency stems from TBB's lightweight task management and the use of static
partitioning for uniform workloads.

### 7.3 TBB vs OpenMP Comparison

The TBB implementation outperforms OpenMP for several reasons:

1. **Lower scheduling overhead:** TBB's `static_partitioner` eliminates the
   runtime scheduling decisions that OpenMP must make, reducing per-iteration
   overhead.

2. **Better cache behavior:** TBB's partitioning strategy may result in more
   cache-friendly memory access patterns for certain workloads.

3. **Efficient thread-local storage:** TBB's internal handling of thread-local
   data reduces false sharing compared to the explicit vector-of-vectors
   approach in the OpenMP version.

4. **Reduced synchronization:** The static partitioner ensures that each
   thread's workload is determined at compile time, eliminating runtime
   synchronization points.

### 7.4 Scalability Analysis

The algorithm demonstrates excellent scalability with problem size. Using linear
approximation, the execution time on the test machine for the TBB version can be
estimated by the formula: `time_tbb (ms) = 0.0000825 × N + 2.45`

Comparing with the sequential formula `time_seq (ms) = 0.0005 × N + 20.25`:

- The TBB implementation achieves approximately 6.1 times lower slope
  coefficient
- Constant overhead is reduced from roughly 20 milliseconds to about 2.5
  milliseconds

Comparing with the OpenMP formula `time_omp (ms) = 0.00011 × N + 3.85`:

- TBB achieves 25% lower slope coefficient
- Constant overhead is reduced by approximately 36%

### 7.5 Strong and Weak Scaling

**Strong scaling** (fixed problem size of 100 million elements, varying
threads):

| Threads | Time (ms) | Speedup | Efficiency |
| ------- | --------- | ------- | ---------- |
| 1       | 52430     | 1.00x   | 100%       |
| 2       | 26890     | 1.95x   | 97.5%      |
| 4       | 13850     | 3.79x   | 94.8%      |
| 8       | 8250      | 6.35x   | 79.4%      |

The implementation achieves near-linear scaling up to 4 threads and maintains
excellent efficiency at 8 threads, demonstrating TBB's ability to effectively
utilize available hardware resources.

## 8. Conclusions

A parallel LSD radix sort for double-precision numbers has been successfully
implemented and validated using Intel Threading Building Blocks. The algorithm
achieves speedups of 6.47x on 8 hardware threads for large datasets (100 million
elements), outperforming the OpenMP implementation by 33% and demonstrating
superior parallel efficiency. The key innovations include:

- A bit transformation technique that eliminates separate handling of negative
  numbers, simplifying the algorithm
- TBB's static partitioning strategy for predictable, low-overhead parallel
  execution
- Efficient histogram construction and scatter phases with thread-local storage
- Careful selection of partitioning strategies based on workload characteristics

The TBB implementation demonstrates that high-level parallel programming
frameworks can achieve better performance than lower-level threading models when
used appropriately, particularly for regular, predictable workloads. The main
limitation remains the O(n) additional memory requirement, though this is
inherent to LSD radix sort implementations.

## 9. References

1. [Сортировки. Из курса "Параллельные численые методы" Сиднев А.А., Сысоев А.В., Мееров И.Б.](http://www.hpcc.unn.ru/file.php?id=458)

2. [Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2009). Introduction to Algorithms (3rd ed.). MIT Press. (Chapter 8: Sorting in Linear Time)](https://ressources.unisciel.fr/algoprog/s00aaroot/aa00module1/res/%5BCormen-AL2011%5DIntroduction_To_Algorithms-A3.pdf)

3. [Knuth, D. E. (1998). The Art of Computer Programming, Volume 3: Sorting and Searching (2nd ed.). Addison-Wesley.](<https://kolegite.com/EE_library/books_and_lectures/Програмиране/The_Art_of_Computer_Programming/The%20Art%20of%20Computer%20Programming%20Volume%203%20Sorting%20and%20Searching%20(Donald%20E.%20Knuth)%20(z-lib.org).pdf>)

4. [Intel Threading Building Blocks Developer Guide, 2021.11.0](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb-documentation.html)
