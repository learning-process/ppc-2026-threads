# Radix sort of `double`s with simple merge (OpenMP parallelization)

- Student: Митяева Дарья Викторовна, group 3823Б1ФИ2
- Technology: OMP
- Variant: 19

## 1. Introduction

Sorting is a fundamental operation in computer science with applications across
all domains of computing. This project implements a parallel Least Significant
Digit (LSD) radix sort specifically designed for double-precision floating-point
numbers using OpenMP for shared-memory parallelization. The algorithm leverages
the byte representation of doubles to achieve linear-time complexity while
utilizing multiple threads to accelerate the sorting process. A key enhancement
in this implementation is the transformation of IEEE 754 double representation
into a sortable unsigned integer format, eliminating the need for separate
handling of negative numbers.

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

## 4. Parallelization Scheme (OpenMP)

The OpenMP implementation parallelizes the radix sort through the following
strategies:

### 4.1 Bit Transformation

Instead of handling negative numbers separately, this implementation uses a
transformation function that maps IEEE 754 doubles to a sortable unsigned
integer representation. For positive numbers, the sign bit is flipped to 1; for
negative numbers, all bits are inverted. This elegant approach ensures that the
integer order matches the floating-point order completely, eliminating the need
for separate negative/positive processing paths and simplifying the algorithm
significantly.

### 4.2 Parallel Counting Pass

Each counting sort pass is parallelized using a three-phase approach:

**Phase 1 – Parallel histogram construction:** Threads process disjoint chunks
of the input array, with each thread building its own local histogram of byte
frequencies (256 buckets per thread). This avoids contention on shared counters.

**Phase 2 – Sequential prefix sum aggregation:** The thread-local histograms are
combined into global prefix sums. Although this phase is sequential, it operates
on only 256 × T elements (where T is the number of threads), which is negligible
compared to the main data processing.

**Phase 3 – Parallel scatter:** Using the computed prefix sums, threads write
elements to their final positions in parallel. Each thread maintains its own
position pointers, ensuring no write conflicts.

### 4.3 Work Distribution

For an array of N elements and T threads, the work is distributed by dividing
the array into T contiguous chunks of approximately N/T elements. The last
thread handles any remainder elements. This contiguous partitioning ensures good
cache locality and minimizes false sharing between threads.

### 4.4 Memory Management

The implementation uses double buffering with two arrays – a source and a
destination – alternating between them after each counting pass. This approach
avoids repeated memory allocations and uses pointer swapping instead of copying.
The number of threads is dynamically configured through the framework's utility
functions.

## 5. Implementation Details

### File Structure

- `common/include/common.hpp` – Type aliases for input, output, and test data
- `common/include/test_generator.hpp` – Random double vector generation
  utilities
- `omp/include/ops_omp.hpp` – Task class interface for framework integration
- `omp/include/sorter_omp.hpp` – Sorting algorithm interface declaration
- `omp/src/ops_omp.cpp` – Framework wrapper implementation (validation,
  preprocessing, execution)
- `omp/src/sorter_omp.cpp` – Core sorting algorithm with OpenMP directives

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
  the three-phase parallel counting process.

- **Sort** – The main sorting routine that orchestrates all passes. It
  transforms the input doubles to sortable integers, performs eight counting
  passes (one per byte), and then converts the sorted integers back to doubles.

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
- **Toolchain:** GCC 14.2.0 x86-64-linux-gnu, build type Release
- **Environment:** OpenMP parallel execution, 8 threads
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

### 7.2 Performance

The following table shows execution times for various input sizes compared to
the sequential baseline:

| Mode | Count       | Time (ms) | Speedup vs Seq |
| ---- | ----------- | --------- | -------------- |
| seq  | 10          | 14        | 1.00x          |
| omp  | 10          | 4         | 3.50x          |
| seq  | 100         | 17        | 1.00x          |
| omp  | 100         | 5         | 3.40x          |
| seq  | 1,000       | 16        | 1.00x          |
| omp  | 1,000       | 4         | 4.00x          |
| seq  | 10,000      | 18        | 1.00x          |
| omp  | 10,000      | 5         | 3.60x          |
| seq  | 100,000     | 77        | 1.00x          |
| omp  | 100,000     | 18        | 4.28x          |
| seq  | 1,000,000   | 495       | 1.00x          |
| omp  | 1,000,000   | 108       | 4.58x          |
| seq  | 10,000,000  | 5138      | 1.00x          |
| omp  | 10,000,000  | 1072      | 4.79x          |
| seq  | 100,000,000 | 53375     | 1.00x          |
| omp  | 100,000,000 | 10984     | 4.86x          |

**Analysis:** The OpenMP parallel implementation demonstrates excellent scaling
with input size, achieving speedups between 3.4x and 4.86x compared to the
sequential baseline. The speedup improves with larger datasets as the parallel
overhead becomes amortized over more work.

Key observations:

- **Small arrays (under 10,000 elements):** Speedup is slightly lower
  (3.4x–4.0x) due to OpenMP thread creation and synchronization overhead
  dominating the execution time.

- **Medium arrays (100,000 – 1,000,000 elements):** Speedup reaches 4.28x–4.58x
  as work distribution becomes more efficient and the overhead becomes
  negligible.

- **Large arrays (10M – 100M elements):** Speedup stabilizes around 4.79x–4.86x,
  approaching the theoretical maximum for 8 hardware threads. The primary
  limiting factor is memory bandwidth, as each pass reads and writes the entire
  dataset.

The parallel efficiency (speedup divided by thread count) ranges from 42.5% to
60.8%, which is excellent for a memory-bound algorithm like radix sort. The
efficiency increases with problem size, reaching its peak at the largest
dataset.

### 7.3 Scalability Analysis

The algorithm demonstrates near-linear scaling with problem size. Using linear
approximation, the execution time on the test machine for the parallel version
can be estimated by the formula: `time_omp (ms) = 0.00011 × N + 3.85`

Comparing with the sequential formula `time_seq (ms) = 0.0005 × N + 20.25`:

- The parallel implementation achieves approximately 4.5 times lower slope
  coefficient
- Constant overhead is reduced from roughly 20 milliseconds to about 4
  milliseconds due to the more efficient bit transformation approach that
  eliminates separate negative number processing

### 7.4 Comparison with Sequential Version

The OpenMP version offers several advantages beyond raw speed:

1. **Simplified negative number handling:** The bit transformation approach
   eliminates the need for separate sorting of negative and positive numbers,
   reducing code complexity and conditional branches.

2. **Better memory locality:** Parallel threads access contiguous memory
   regions, improving cache utilization and reducing cache misses.

3. **Reduced constant factors:** The transformation approach uses a more uniform
   processing path with less conditional logic, contributing to the lower
   constant overhead observed in measurements.

## 8. Conclusions

A parallel LSD radix sort for double-precision numbers has been successfully
implemented and validated using OpenMP. The algorithm achieves speedups of 4.86x
on 8 hardware threads for large datasets (100 million elements), demonstrating
excellent parallel efficiency. The key innovations include:

- A bit transformation technique that eliminates separate handling of negative
  numbers, simplifying the algorithm
- Parallel histogram construction using thread-local counters to avoid
  contention
- An efficient scatter phase using thread-local position pointers for
  conflict-free writes

The implementation serves as a strong baseline for further parallelization using
MPI for distributed memory systems or hybrid approaches combining OpenMP with
MPI. The main limitation remains the O(n) additional memory requirement, though
this is inherent to LSD radix sort implementations and could be addressed in
future work through in-place radix sort techniques.

## 9. References

1. [Сортировки. Из курса "Параллельные численые методы" Сиднев А.А., Сысоев А.В., Мееров И.Б.](http://www.hpcc.unn.ru/file.php?id=458)

2. [Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2009). Introduction to Algorithms (3rd ed.). MIT Press. (Chapter 8: Sorting in Linear Time)](https://ressources.unisciel.fr/algoprog/s00aaroot/aa00module1/res/%5BCormen-AL2011%5DIntroduction_To_Algorithms-A3.pdf)

3. [Knuth, D. E. (1998). The Art of Computer Programming, Volume 3: Sorting and Searching (2nd ed.). Addison-Wesley.](<https://kolegite.com/EE_library/books_and_lectures/Програмиране/The_Art_of_Computer_Programming/The%20Art%20of%20Computer%20Programming%20Volume%203%20Sorting%20and%20Searching%20(Donald%20E.%20Knuth)%20(z-lib.org).pdf>)

4. [OpenMP Application Programming Interface Specification Version 5.0](https://www.openmp.org/spec-html/5.0/openmp50.html)
