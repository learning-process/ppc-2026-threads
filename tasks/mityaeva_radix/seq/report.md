# Radix sort of `double`s with simple merge

- Student: Митяева Дарья Викторовна, group 3823Б1ФИ2
- Technology: SEQ
- Variant: 19

## 1. Introduction

Sorting is a fundamental operation in computer science with applications across all domains of computing. This project
implements a sequential Least Significant Digit (LSD) radix sort specifically designed for double-precision
floating-point numbers. The algorithm leverages the byte representation of doubles to achieve linear-time complexity,
providing an efficient alternative to comparison-based sorting algorithms for large datasets. A key enhancement in this
implementation is the separate handling of negative and positive numbers to address the IEEE 754 representation quirk
where negative doubles have a different byte ordering.

## 2. Problem Statement

- **Input:** A vector of double-precision floating-point numbers of arbitrary length `N`.
- **Output:** The same vector sorted in non-decreasing order.

**Constraints:** The input vector must not be empty. The algorithm must handle all possible double values including
positive/negative zero.

## 3. Baseline Algorithm (Sequential)

The LSD radix sort processes numbers digit by digit from least significant to most significant. For double values (8
bytes = 64 bits), the algorithm:

1. Interprets each double as an array of 8 unsigned bytes.
2. Performs counting sort on each byte position (0 to 7).
3. Alternates between original and auxiliary arrays to avoid unnecessary copying.
4. Maintains stability throughout all passes to ensure correct final ordering.

The counting sort for each byte:

- Builds a histogram of byte values (256 buckets).
- Computes prefix sums to determine final positions.
- Places elements in their sorted positions maintaining stability.

Due to the `IEEE 754` representation, negative numbers require special handling: their byte representation is inverted
to maintain correct numerical order.

## 4. Parallelization Scheme

This is a sequential implementation, designed as a baseline for future parallel comparisons. The algorithm serves as a
reference point for evaluating the performance gains of parallel versions using MPI, OpenMP, TBB, or STL.

## 5. Implementation Details

### File Structure

- `common/include/common.hpp` - Type aliases (InType, OutType, TestType)
- `common/include/test_generator.hpp` - Random double vector generation
- `seq/include/ops_seq.hpp` - Task class interface for framework integration
- `seq/include/sorter_seq.hpp` - Sorting algorithm interface
- `seq/src/ops_seq.cpp` - Framework wrapper implementation
- `seq/src/sorter_seq.cpp` - Core sorting algorithm

### Key Functions

- `SorterSeq::CountingSortAsc` - Performs counting sort in ascending order for a specific byte.
- `SorterSeq::CountingSortDesc` - Performs counting sort in descending order for a specific byte (used for negative
  numbers).
- `SorterSeq::LSDSortDouble` - Separates negative and positive numbers, applies appropriate sorting to each group, then
  merges them.

### Negative Numbers Handling

The implementation first separates negative numbers (which are sorted in descending order to account for their inverted
bit representation) from non-negative numbers (sorted in ascending order). After sorting both groups independently, they
are merged with negative numbers placed before positives, maintaining the correct overall order.

## 6. Experimental Setup

- **Hardware/OS:** Intel Core i7-1165G7 @ 2.80GHz (4 cores, 8 threads), 16GB RAM, Ubuntu 22.04 via WSL2 under Windows 10
  (build 2H22)
- **Toolchain:** GCC 14.2.0 x86-64-linux-gnu, build type Release
- **Environment:** Sequential execution, single thread
- **Data:** Random doubles uniformly distributed between -0.5 and 0.5, generated with fixed seed for reproducibility

## 7. Results and Discussion

### 7.1 Correctness

Correctness was verified through:

- Comparison with std::ranges::is_sorted results for multiple random datasets.
- Edge case testing: single element, duplicate values, already sorted arrays, reverse sorted arrays.
- Special value handling: positive and negative zero.
- Extensive testing with mixed positive and negative numbers to ensure proper ordering across the sign boundary.

### 7.2 Performance

The following table shows execution times for various input sizes:

| Mode | Count       | Time (ms) |
| ---- | ----------- | --------- |
| seq  | 10          | 14        |
| seq  | 100         | 17        |
| seq  | 1,000       | 16        |
| seq  | 10,000      | 18        |
| seq  | 100,000     | 77        |
| seq  | 1,000,000   | 495       |
| seq  | 10,000,000  | 5138      |
| seq  | 100,000,000 | 53375     |

**Analysis:** The algorithm demonstrates excellent linear scaling with input size. The linear correlation coefficient
between input size and execution time is `0.9996`, confirming the theoretical `O(n)` complexity of radix sort. Using
linear approximation, the execution time on the test machine can be estimated by the formula:
`time (ms) = 0.0005 × N + 20.2538`, where `N` is the number of elements to sort. This formula is obtained through
approximation and should be considered approximate, as actual performance may vary based on data distribution, memory
hierarchy effects, and system load.

Observations:

- The overhead for small arrays (under 10,000 elements) is relatively constant at around 14-18 ms, dominated by function
  call overhead and vector allocations.
- Performance scales linearly once arrays exceed 100,000 elements.
- The algorithm handles 100 million doubles (800 MB of data) in under 1 minute, demonstrating excellent efficiency.
- The separate handling of negative numbers adds minimal overhead while ensuring correctness.
- Memory bandwidth becomes a noticeable bottleneck for the largest dataset, though the impact is less pronounced than in
  comparison-based sorts.

## 8. Conclusions

A sequential LSD radix sort for double-precision numbers has been successfully implemented and validated, with proper
handling of negative numbers. The algorithm achieves linear time complexity with a correlation coefficient of 0.9996,
making it highly efficient for large-scale sorting tasks. The implementation successfully addresses the `IEEE 754`
representation challenge for negative doubles through separate sorting paths. It serves as a solid baseline for future
parallel implementations using various parallel programming technologies. The main limitation is the O(n) additional
memory requirement, which could be addressed in future work through in-place radix sort techniques.

## 9. References

1. [Сортировки. Из курса "Параллельные численые методы" Сиднев А.А., Сысоев А.В., Мееров И.Б.](http://www.hpcc.unn.ru/file.php?id=458)
2. [Cormen, T. H., Leiserson, C. E., Rivest, R. L., & Stein, C. (2009). Introduction to Algorithms (3rd ed.). MIT Press. (Chapter 8: Sorting in Linear Time)](https://ressources.unisciel.fr/algoprog/s00aaroot/aa00module1/res/%5BCormen-AL2011%5DIntroduction_To_Algorithms-A3.pdf)
3. [Knuth, D. E. (1998). The Art of Computer Programming, Volume 3: Sorting and Searching (2nd ed.). Addison-Wesley.](<https://kolegite.com/EE_library/books_and_lectures/Програмиране/The_Art_of_Computer_Programming/The%20Art%20of%20Computer%20Programming%20Volume%203%20Sorting%20and%20Searching%20(Donald%20E.%20Knuth)%20(z-lib.org).pdf>)
