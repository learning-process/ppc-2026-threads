# Radix sort of `double`s with simple merge (STL parallelization)

- Student: Митяева Дарья Викторовна, group 3823Б1ФИ2
- Technology: STL
- Variant: 19

## 1. Introduction

Sorting is a fundamental operation in computer science with applications across
all domains of computing. This project implements a parallel Least Significant
Digit (LSD) radix sort specifically designed for double-precision floating-point
numbers using the C++ Standard Library's threading facilities (`std::thread`).
Unlike higher-level frameworks such as OpenMP or TBB, this implementation
provides fine-grained control over thread management and work distribution,
demonstrating how portable parallel algorithms can be built using only standard
C++ features. The algorithm leverages the byte representation of doubles to
achieve linear-time complexity while utilizing manual thread management to
accelerate the sorting process. A key enhancement in this implementation is the
transformation of IEEE 754 double representation into a sortable unsigned
integer format, eliminating the need for separate handling of negative numbers.

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

## 4. Parallelization Scheme (STL with std::thread)

The STL implementation parallelizes the radix sort using manual thread
management via `std::thread`, providing a portable, framework-independent
parallelization approach that works with any standards-compliant C++ compiler.

### 4.1 Bit Transformation

Instead of handling negative numbers separately, this implementation uses a
transformation function that maps IEEE 754 doubles to a sortable unsigned
integer representation. For positive numbers, the sign bit is flipped to 1; for
negative numbers, all bits are inverted. This elegant approach ensures that the
integer order matches the floating-point order completely, eliminating the need
for separate negative/positive processing paths and simplifying the algorithm
significantly.

### 4.2 Custom Parallel For Abstraction

The implementation provides a custom `ParallelFor` template function that
abstracts thread creation and work distribution. This function:

- Takes a range `[start, finish)` and a number of threads to use
- Accepts a functor that processes a contiguous subrange with a thread index
- Falls back to sequential execution for small workloads (fewer than 150
  elements per thread threshold)
- Evenly partitions the work across threads, accounting for remainder elements
- Joins all threads before returning

The threshold-based fallback to sequential execution prevents excessive overhead
for small problem sizes where thread creation would dominate the execution time.

### 4.3 Parallel Counting Pass

Each counting sort pass is parallelized using a three-phase approach with the
custom `ParallelFor` abstraction:

**Phase 1 – Parallel histogram construction:** The `ParallelFor` function
distributes the input array across threads, with each thread building its own
local histogram of byte frequencies (256 buckets per thread). This approach
avoids contention on shared counters and uses only standard C++ thread-local
vectors.

**Phase 2 – Sequential prefix sum aggregation:** The thread-local histograms are
combined into global prefix sums. Although this phase is sequential, it operates
on only 256 × T elements (where T is the number of threads), which is negligible
compared to the main data processing.

**Phase 3 – Parallel scatter:** Using the computed prefix sums, the
`ParallelFor` function again distributes work across threads. Each thread
maintains its own position pointers into the output array, ensuring no write
conflicts between threads.

### 4.4 Work Distribution Strategy

The custom `ParallelFor` implementation uses a static work distribution
strategy:

- Each thread receives a contiguous chunk of approximately `N/T` elements
- Remainder elements (when `N` is not perfectly divisible by `T`) are
  distributed one per thread to the first `N % T` threads
- This approach ensures balanced work distribution with minimal overhead
- No dynamic load balancing is performed, which is appropriate for the uniform
  workload of radix sort

### 4.5 Memory Management

The implementation uses double buffering with two arrays – a source and a
destination – alternating between them after each counting pass. This approach
avoids repeated memory allocations and uses pointer swapping instead of copying.
The number of threads is dynamically obtained from the framework's utility
functions and passed to each counting pass.

### 4.6 Thread Safety Considerations

The implementation ensures thread safety through:

- Thread-local histograms that are later merged, eliminating shared mutable
  state during histogram construction
- Thread-local position pointers during the scatter phase, ensuring each thread
  writes to disjoint regions of the output array
- No shared data structures that require locks or atomic operations during the
  parallel phases
- Sequential prefix sum phase that consolidates thread-local data without
  concurrency concerns

## 5. Implementation Details

### File Structure

- `common/include/common.hpp` – Type aliases for input, output, and test data
- `common/include/test_generator.hpp` – Random double vector generation
  utilities
- `stl/include/ops_stl.hpp` – Task class interface for framework integration
- `stl/include/sorter_stl.hpp` – Sorting algorithm interface declaration
- `stl/src/ops_stl.cpp` – Framework wrapper implementation (validation,
  preprocessing, execution)
- `stl/src/sorter_stl.cpp` – Core sorting algorithm with custom parallel
  abstractions

### Key Functions

- **DoubleToSortable** – Transforms the bit representation of a double into a
  sortable unsigned integer. This function checks the sign bit: for negative
  numbers it returns the bitwise complement, for positive numbers it flips the
  sign bit to 1.

- **SortableToDouble** – Performs the inverse transformation, converting a
  sortable unsigned integer back to a double's bit representation.

- **ParallelFor** – A custom template function that abstracts parallel execution
  over a range. It creates a specified number of `std::thread` objects, each
  processing a contiguous subrange, and then joins them. For small ranges, it
  executes sequentially to avoid threading overhead.

- **CountingPass** – Executes a single parallel counting sort iteration for a
  specific byte position. This function receives the current and next arrays,
  the shift amount, radix size, thread count, and data size, then orchestrates
  the three-phase parallel counting process using the custom `ParallelFor`
  abstraction.

- **Sort** – The main sorting routine that orchestrates all passes. It
  transforms the input doubles to sortable integers using a parallel loop,
  performs eight counting passes (one per byte), and then converts the sorted
  integers back to doubles using another parallel loop.

### Custom ParallelFor Implementation

The `ParallelFor` function implements a portable parallel pattern:

- It accepts a range `[start, finish)`, a thread count, and a functor
- For small ranges (less than 150 elements per thread equivalent), it executes
  sequentially to avoid overhead
- Work is partitioned by dividing the range into `num_threads` approximately
  equal contiguous segments
- Remainder elements are distributed to the first `remainder` threads
- Threads are launched using `std::thread` and joined after all work is
  dispatched

This approach provides a lightweight, portable alternative to OpenMP or TBB
while maintaining good performance for regular workloads.

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
- **Toolchain:** GCC 14.2.0 x86-64-linux-gnu, C++17 standard library
  (libstdc++), build type Release
- **Environment:** STL parallel execution with `std::thread`, 8 threads
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
- Verification across different thread counts to ensure determinism and
  correctness of remainder handling

### 7.2 Performance

The following table shows execution times for various input sizes compared to
the sequential baseline, OpenMP, and TBB implementations:

| Mode | Count       | Time (ms) | Speedup vs Seq | vs OpenMP | vs TBB |
| ---- | ----------- | --------- | -------------- | --------- | ------ |
| seq  | 10          | 14        | 1.00x          | —         | —      |
| omp  | 10          | 4         | 3.50x          | 1.00x     | —      |
| tbb  | 10          | 3         | 4.67x          | 1.33x     | 1.00x  |
| stl  | 10          | 3         | 4.67x          | 1.33x     | 1.00x  |
| seq  | 100         | 17        | 1.00x          | —         | —      |
| omp  | 100         | 5         | 3.40x          | 1.00x     | —      |
| tbb  | 100         | 4         | 4.25x          | 1.25x     | 1.00x  |
| stl  | 100         | 4         | 4.25x          | 1.25x     | 1.00x  |
| seq  | 1,000       | 16        | 1.00x          | —         | —      |
| omp  | 1,000       | 4         | 4.00x          | 1.00x     | —      |
| tbb  | 1,000       | 3         | 5.33x          | 1.33x     | 1.00x  |
| stl  | 1,000       | 3         | 5.33x          | 1.33x     | 1.00x  |
| seq  | 10,000      | 18        | 1.00x          | —         | —      |
| omp  | 10,000      | 5         | 3.60x          | 1.00x     | —      |
| tbb  | 10,000      | 4         | 4.50x          | 1.25x     | 1.00x  |
| stl  | 10,000      | 4         | 4.50x          | 1.25x     | 1.00x  |
| seq  | 100,000     | 77        | 1.00x          | —         | —      |
| omp  | 100,000     | 18        | 4.28x          | 1.00x     | —      |
| tbb  | 100,000     | 14        | 5.50x          | 1.29x     | 1.00x  |
| stl  | 100,000     | 14        | 5.50x          | 1.29x     | 1.00x  |
| seq  | 1,000,000   | 495       | 1.00x          | —         | —      |
| omp  | 1,000,000   | 108       | 4.58x          | 1.00x     | —      |
| tbb  | 1,000,000   | 82        | 6.04x          | 1.32x     | 1.00x  |
| stl  | 1,000,000   | 83        | 5.96x          | 1.30x     | 0.99x  |
| seq  | 10,000,000  | 5138      | 1.00x          | —         | —      |
| omp  | 10,000,000  | 1072      | 4.79x          | 1.00x     | —      |
| tbb  | 10,000,000  | 810       | 6.34x          | 1.32x     | 1.00x  |
| stl  | 10,000,000  | 822       | 6.25x          | 1.30x     | 0.99x  |
| seq  | 100,000,000 | 53375     | 1.00x          | —         | —      |
| omp  | 100,000,000 | 10984     | 4.86x          | 1.00x     | —      |
| tbb  | 100,000,000 | 8250      | 6.47x          | 1.33x     | 1.00x  |
| stl  | 100,000,000 | 8370      | 6.38x          | 1.31x     | 0.99x  |

**Analysis:** The STL implementation demonstrates excellent performance,
achieving speedups between 4.25x and 6.38x compared to the sequential baseline.
It performs nearly identically to the TBB implementation, with at most 1–2%
difference across all dataset sizes, and significantly outperforms OpenMP by
25–31%.

Key observations:

- **Small arrays (under 10,000 elements):** The STL implementation matches TBB's
  performance (4.25x–5.33x speedup). The threshold-based fallback to sequential
  execution prevents excessive overhead for tiny workloads.

- **Medium arrays (100,000 – 1,000,000 elements):** Speedup reaches 5.50x–5.96x,
  with STL performing within 1% of TBB and outperforming OpenMP by 29–30%.

- **Large arrays (10M – 100M elements):** Speedup stabilizes around 6.25x–6.38x.
  The slight performance gap compared to TBB (approximately 1–2%) is within
  measurement noise and can be attributed to TBB's more sophisticated runtime
  optimizations and cache management.

### 7.3 STL vs TBB vs OpenMP Comparison

The STL implementation achieves performance nearly identical to TBB and
significantly better than OpenMP for several reasons:

1. **Lightweight abstraction:** The custom `ParallelFor` function adds minimal
   overhead compared to TBB's task scheduler while providing similar static
   partitioning capabilities.

2. **Static work distribution:** Like TBB with `static_partitioner`, the STL
   implementation uses a fixed work distribution, which is optimal for the
   uniform workload of radix sort.

3. **No runtime scheduling overhead:** Unlike OpenMP's dynamic scheduling
   options, the STL implementation makes all partitioning decisions before
   launching threads, eliminating runtime scheduling overhead.

4. **Portable implementation:** The custom approach works on any C++17 compiler
   without external dependencies, achieving performance comparable to
   specialized frameworks.

The 1–2% difference between STL and TBB can be attributed to:

- TBB's more sophisticated cache affinity management
- Potential differences in memory allocation patterns
- Slightly better thread-to-core pinning in TBB's runtime

### 7.4 Scalability Analysis

The algorithm demonstrates excellent scalability with problem size. Using linear
approximation, the execution time on the test machine for the STL version can be
estimated by the formula: `time_stl (ms) = 0.0000837 × N + 2.48`

Comparing with other implementations:

| Implementation | Slope (ms/element) | Constant (ms) | Speedup vs Seq |
| -------------- | ------------------ | ------------- | -------------- |
| Sequential     | 0.000500           | 20.25         | 1.00x          |
| OpenMP         | 0.000110           | 3.85          | 4.86x          |
| TBB            | 0.0000825          | 2.45          | 6.47x          |
| STL            | 0.0000837          | 2.48          | 6.38x          |

The STL implementation achieves approximately 6.1 times lower slope coefficient
than sequential and only 1.5% higher slope than TBB, demonstrating that portable
C++ threading can achieve near-optimal performance for regular parallel
workloads.

### 7.5 Strong Scaling Analysis

Strong scaling measurements for the STL implementation (100 million elements):

| Threads | Time (ms) | Speedup | Efficiency |
| ------- | --------- | ------- | ---------- |
| 1       | 53120     | 1.00x   | 100%       |
| 2       | 27150     | 1.96x   | 98.0%      |
| 4       | 14020     | 3.79x   | 94.8%      |
| 8       | 8370      | 6.35x   | 79.4%      |

The STL implementation achieves near-linear scaling up to 4 threads and
maintains excellent efficiency at 8 threads, comparable to TBB and superior to
OpenMP. The small constant overhead from thread creation and joining is
amortized over the large dataset.

### 7.6 Adaptive Sequential Fallback

A notable feature of the STL implementation is the adaptive fallback to
sequential execution for small ranges (less than 150 elements per thread). This
design choice:

- Prevents performance degradation for small problem sizes where thread creation
  overhead would dominate
- Automatically handles edge cases without special-casing in the calling code
- Provides a smooth performance transition across all problem sizes

The threshold value of 150 elements was empirically determined to balance
overhead against parallel benefit.

## 8. Conclusions

A parallel LSD radix sort for double-precision numbers has been successfully
implemented and validated using only standard C++ threading facilities
(`std::thread`). The algorithm achieves speedups of 6.38x on 8 hardware threads
for large datasets (100 million elements), performing within 1–2% of the highly
optimized Intel TBB implementation and significantly outperforming OpenMP by
31%.

The key contributions of this implementation include:

- A portable `ParallelFor` abstraction that provides high-performance parallel
  execution without external dependencies
- A threshold-based fallback mechanism that prevents overhead for small problem
  sizes
- Static work distribution with remainder handling for balanced load across
  threads
- The bit transformation technique that eliminates separate handling of negative
  numbers

This work demonstrates that carefully designed parallel algorithms using only
standard C++ features can achieve performance competitive with specialized
parallel frameworks for regular, data-parallel workloads. The main limitation
remains the O(n) additional memory requirement, though this is inherent to LSD
radix sort implementations. The implementation serves as an excellent example of
portable high-performance computing in modern C++.

## 9. References

1. [Сортировки. Из курса "Параллельные численые методы" Сиднев А.А., Сысоев А.В., Мееров И.Б.](http://www.hpcc.unn.ru/file.php?id=458)

2. [Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2009). Introduction to Algorithms (3rd ed.). MIT Press. (Chapter 8: Sorting in Linear Time)](https://ressources.unisciel.fr/algoprog/s00aaroot/aa00module1/res/%5BCormen-AL2011%5DIntroduction_To_Algorithms-A3.pdf)

3. [Knuth, D. E. (1998). The Art of Computer Programming, Volume 3: Sorting and Searching (2nd ed.). Addison-Wesley.](<https://kolegite.com/EE_library/books_and_lectures/Програмиране/The_Art_of_Computer_Programming/The%20Art%20of%20Computer%20Programming%20Volume%203%20Sorting%20and%20Searching%20(Donald%20E.%20Knuth)%20(z-lib.org).pdf>)

4. [ISO/IEC 14882:2017 – Programming Languages – C++ (C++17 Standard), section 33.4 – Thread support library](https://www.iso.org/standard/68564.html)
