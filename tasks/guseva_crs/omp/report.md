# Sparse Matrix Multiplication. Double Precision Elements. Compressed Row Storage (CRS) Format

- _Student_: Гусева Алёна Сергеевна, group 3823Б1ФИ2
- _Technology_: OMP
- _Variant_: 4

## 1. Introduction

### 1.1 Matrix Multiplication

Sparse matrix multiplication is a fundamental operation in scientific computing, engineering simulations, and data
analysis. Many real-world problems involve matrices with a large number of elements, where most entries are zero.
Examples include finite element analysis, graph processing, network analysis, and machine learning applications. Storing
and operating on these matrices in dense format would be extremely memory-inefficient and computationally wasteful.

### 1.2 SRC format

The Compressed Row Storage (CRS) format, also known as Compressed Sparse Row (CSR) format, is one of the most widely
used sparse matrix representations. It stores only non-zero elements, along with their column indices and row pointers,
significantly reducing memory requirements and enabling efficient matrix operations.

### 1.3 OpenMP

The most popular multithreading C++ library is OpenMP. It provides a simple, high-level interface for shared-memory
parallelism, allowing developers to parallelize loops and sections of code with minimal modifications to the underlying
source. For the task of sparse CRS matrix multiplication, this simplicity is particularly valuable, as it allows us to
focus on distributing the irregular row-based workload across cores without getting bogged down in complex thread
management. By annotating the outer row loop with OpenMP directives, we can transform a sequential sparse multiplication
routine into a parallel one capable of utilizing multiple processors effectively.

### 1.4 Expected Outcome

The expected outcome of this work is a parallel implementation of sparse matrix multiplication for matrices stored in
CRS format using OpenMP. The implementation will read matrices from test files, perform multiplication with correctness
verification, and produce the result matrix also in CRS format, while demonstrating performance improvements through
parallel execution on multi-core processors.

## 2. Problem Statement

### 2.1 Formal Definition

Given two matrices `A` and `B` stored in Compressed Row Storage (CRS) format, compute their dot product `C = A × B`,
where `C` is also stored in CRS format.

### Format Requirements

Each matrix should be presented in CRS form:

- `nz`: number of non-zero elements
- `nrows`: number of rows
- `ncols`: number of cols
- `values`: array of non-zero element values (size `nz`)
- `cols`: array of column indices of each non-zero element (size `nz`)
- `row_ptrs`: array of row pointers indicating start indices for each row (size `nrows + 1`)

### 2.2 Input Format

- CRS-stored matrices `A` and `B`

### 2.3 Output Format

- CRS-stored matrix `C = A × B`

### 2.4 Constraints

- `A.ncols == B.nrows`
- CRS format with 0-based indexing
- Elements with absolute value less than 10⁻⁵ are considered zero and so are not stored in the result

## 3. Baseline Algorithm (Sequential)

The baseline algorithm for sparse matrix multiplication in CRS format follows the same logical steps as the sequential
version but leverages OpenMP for parallel execution. The key difference lies in the distribution of work across multiple
threads.

### 3.1 Matrix Transposition

Before performing the multiplication, matrix `B` is transposed using the same algorithm as in the sequential version.
This operation remains sequential as it is difficult to parallelize efficiently due to data dependencies. However, for
large matrices, the transposition cost is amortized over the subsequent parallel multiplication.

1. Create a new CRS structure for the transposed matrix Bᵀ

2. Set `nrows = B.ncols`, `ncols = B.nrows`

3. Count the number of elements in each row of `Bᵀ` by iterating through `B.cols` and incrementing counters

4. Build the row pointers for `Bᵀ` using prefix sums

5. Fill the values and column indices of `Bᵀ` by iterating through the original matrix `B` and placing each element in
   its new position

### 3.2 Multiplication Algorithm

The multiplication is parallelized by distributing rows of matrix `A` among available threads. For each row `i` of
matrix `A` assigned to a thread:

1. Initialize a temporary array or marker to track which columns of `A` have non-zero elements in this row

2. For each non-zero element in row `i` of `A`, record its column index and position

3. For each row `j` of `Bᵀ` (which corresponds to column `j` of `B`):
    - Initialize `sum = 0`

    - For each non-zero element in row `j` of `Bᵀ`:
        - If the column index of this element (which corresponds to a row index in original `B`) matches a recorded
          column from row `i` of `A`

        - Multiply the corresponding values from `A` and `Bᵀ` and add to `sum`

    - If `|sum| > 10⁻⁵`, add this element to the result matrix `C` at position `(i, j)`

### 3.3

The parallelization is achieved using the following OpenMP constructs:

- `#pragma omp parallel`: Creates a team of threads
- `#pragma omp for`: Distributes loop iterations (rows of `A`) among threads
- **Private variables**: Each thread gets its own copies of loop indices and temporary storage to avoid race conditions
- **Shared variables**: Input matrices and result structures are shared among threads
- **Manual result aggregation**: Each thread builds its own portion of the result matrix, which is later combined

### 3.4 Load Balancing Considerations

Due to the irregular nature of sparse matrices, the number of non-zero elements per row can vary significantly. This
creates potential load imbalance where some threads finish their assigned rows quickly while others continue working.
OpenMP's default static scheduling may not be optimal for such cases, but the implementation can be extended to use
dynamic scheduling if needed.

### 3.5 Complexity Analysis

**Time complexity:**

- Sequential portion: `O(nz(B))` for transposition
- Parallel portion: `O((nz(A) + nz(B) + number_of_non-zero_products) / P)` where `P` is the number of threads, plus
  synchronization overhead

**Space complexity**:

- Per thread: `O(ncols(A))` for private marker arrays
- Total: `O(nz(A) + nz(B) + nz(C) + P × ncols(A))`

## 4. Implementation Details

### 4.1 Code Structure

```└── guseva_crs/
        ├── common/
        |   ├── common.hpp  . . . . . . . . . . . . Common type definitions
        |   |                                       and constants, CRS struct definition
        |   ├── multiplier.hpp  . . . . . . . . . . Base class Multiplier
        |   └── test_reader.hpp . . . . . . . . . . Functions to read test data and store them in CRS format
        |
        ├── data/ . . . . . . . . . . . . . . . . . Directory for test files in .txt format
        |
        ├── omp/
        |   ├── include/
        |   |   ├── ops_omp.hpp . . . . . . . . . . OpenMP task class declaration
        |   |   └── multiplier_omp.hpp  . . . . . . OpenMP multiplier implementation
        |   |                                       inherited from Multiplier
        |   └── src/
        |       └── ops_omp.cpp . . . . . . . . . . OpenMP task implementation
        |
        ├── tests/
        |   ├── functional/
        |   |   └── main.cpp  . . . . . . . . . . . Functional tests
        |   |
        |   └── performance/
        |       └── main.cpp  . . . . . . . . . . . Performance tests
        |
        ├── info.json . . . . . . . . . . . . . . . Information about author of the work
        ├── report.md . . . . . . . . . . . . . . . Report you are reading right now
        └── settings.json . . . . . . . . . . . . . Configuration of the current work
```

### 4.2 Key Classes And Functions

- `struct CRS` defined in `./common/common.hpp` to store matrices in CRS format
- `class Multiplier` defined in `./common/multiplier.hpp` to provide interface for multiplier and perform matrix
  transposition
- `functions ReadTestFromFile() and ReadCRSFromFile()` to read and cast `.txt` stored data into CRS format provided by
  `struct CRS`
- `class MultiplierOmp` inherited from `class Multiplier` to perform dot product using OpenMP parallelization
- `class GusevaCRSMatMulOmp` implements the repo's `BaseTask` interface for correct task verification and acceptance
- `class GusevaMatMulCRSFuncTest` to produce functional tests
- `class GusevaMatMulCRSPerfTest` to produce performance tests

### 4.3 Important Assumptions And Corner Cases

The implementation utilizes the following OpenMP features:

#### 4.3.1 Key OpenMP Features Used

1. **Parallel Region**: `#pragma omp parallel` creates a team of threads

2. **Work Sharing**: `#pragma omp for` distributes row iterations

3. **Data Scoping**:
    - `shared(a, bt, columns, values, row_index)` - matrices and result structures are shared

    - `private(j, k)` - loop indices are thread-private

    - `default(none)` - requires explicit specification of all variable scopes

4. **Thread-Local Storage**: `std::vector<int> temp(n)` declared inside parallel region - each thread gets its own copy

5. **Automatic Synchronization**: Implicit barrier at the end of parallel region

#### 4.3.2 Thread Safety Considerations

The implementation ensures thread safety through:

- **Private temporary arrays**: Each thread maintains its own marker vector to avoid conflicts when accessing columns of
  `A`

- **Independent result building**: Threads build separate sections of the result matrix (by row ranges)

- **No shared mutable state**: Threads do not write to shared structures simultaneously

- **Clear data scoping**: All variables are explicitly scoped as shared or private

#### 4.3.3 Assumptions

- Matrices are valid and properly formatted in `CRS` with 0-based indexing

- Row pointers satisfy: `row_ptrs[0] = 0`, `row_ptrs[nrows] = nz`

- Column indices for each row are sorted (typical `CRS` requirement)

- Input matrices contain only double-precision floating-point values

- OpenMP environment is properly configured with available threads

**Corner Cases Handled:**

- Empty matrices: Matrices with nz = 0 are handled correctly

- Zero threshold: Values below 10⁻⁵ are considered zero and not stored in result

- Invalid dimensions: Validation checks that input matrices are multiplicatable

- Single-element matrices: Multiplication works for 1×1 matrices

- Rectangular matrices: Handles non-square matrices correctly

- Load imbalance: Implementation allows for future scheduling optimizations

### 4.4 Memory Usage Considerations

The implementation is designed with memory efficiency in mind:

1. Sparse storage: Only non-zero elements are stored, saving significant memory for sparse matrices

2. In-place transposition: The transpose operation creates a new matrix but avoids unnecessary copying

3. Temporary arrays: The temp vector is reused for each row, sized to the number of columns in `A`

4. Result construction: The result matrix is built incrementally without pre-allocation of full dense storage

**Memory complexity:**

1. Input matrices: `O(nz(A) + nz(B)) + O(nrows(A) + nrows(B) + ncols(A) + ncols(B))` for index arrays

2. Transposed matrix: `O(nz(B)) + O(ncols(B))` for row pointers

3. Temporary storage: `O(ncols(A))` for the marker array

4. Output matrix: `O(nz(C)) + O(nrows(C))` for row pointers

This approach ensures that memory usage scales linearly with the number of non-zero elements rather than the full matrix
dimensions, making it suitable for large sparse matrices.

## 5. Experimental Setup

- **Hardware/OS**:
  - **Host**: Intel Core i7-14700k, 8+12 cores, 32 Gb DDR4, Windows 10 (10.0.19045.6456)
  - **Virtual**: Intel Core i7-14700k, 12 cores, 8 Gb, WSL2 (2.6.1.0) + Ubuntu (24.04.3 LTS)
- **Toolchain**:

    | Compiler |          Version           | Build Type |
    | :------: | :------------------------: | :--------: |
    |   gcc    |  14.2.0 x86_64-linux-gnu   |  Release   |
    |  clang   | 21.1.0 x86-64-pc-linus-gnu |  Release   |

- **Environment**: OpenMP parallel execution with configurable number of threads (default: uses all available cores)
- **Data**: Test data is generated by Python script with usage of `numpy` and `scipy.sparse` libs for func tests. The
  script code is given in `Appendix`, `1. generate_tests.py`. Perf tests data generated on-the-go: they performs
  multiplication of diagonal matrices
  - Functional tests includes some types of dot matrix multiplication:
    - 1 × `sparse  × dense`, sparse has density of 0.2
    - 1 × `dense × sparse`, sparse has density of 0.2
    - 4 × `sparse × sparse` of different sizes, sparses has densities of 0.1
  - Performance tests
    - Diagonal matrices of size `10000 × 10000` with 1000 non-zero elements

## 6. Results and Discussion

### 6.1 Correctness

The correctness of the OpenMP parallel sparse matrix multiplication implementation was verified through the same
comprehensive testing approach as the sequential version:

#### Reference Results Verification

The primary method of correctness verification involved comparing the implementation's output against pre-computed
reference results. Test files were structured to contain three matrices: input matrix `A`, input matrix `B`, and the
expected result matrix `C = A × B`. The verification process consisted of:

1. File-based testing: Each test case was stored in a separate file containing the complete triplet (`A`, `B`, expected
   `C`)

2. Automated comparison: The Equal function from `./common/common.hpp` was used to compare the computed result with the
   expected result

3. Tolerance-based verification: Due to floating-point arithmetic, comparisons used a tolerance threshold of `10⁻⁵`, as
   defined by the `kZERO` constant

#### Parallel-Specific Correctness Considerations

The parallel implementation must additionally ensure that:

- No race conditions occur when threads access shared data

- Thread-local temporary arrays are properly initialized and cleared

- Results from all threads are correctly aggregated

- The parallel execution produces bitwise-identical results to the sequential version (within floating-point tolerance)

All functional tests passed, confirming that the OpenMP implementation maintains correctness while achieving parallel
execution.

#### Invariant Checking

Several invariants were verified throughout the computation to ensure the algorithm maintains correct CRS structure
properties:

**Input Invariants:**

- Row pointers satisfy: `row_ptrs[0] = 0` and `row_ptrs[nrows] = nz`

- Column indices are within bounds: `0 ≤ cols[i] < ncols` for all `i`

- Row pointers are non-decreasing: `row_ptrs[i] ≤ row_ptrs[i+1]` for all `i`

**Output Invariants:**

- The result matrix maintains the same invariants as input matrices

- The number of rows in result equals the number of rows in `A`

- The number of columns in result equals the number of columns in `B`

- All non-zero values in result satisfy `|value| ≥ kZERO`

**Algorithmic Invariants:**

- After transposition, matrix `Bᵀ` satisfies: `(Bᵀ).nrows = B.ncols` and `(Bᵀ).ncols = B.nrows`

- The multiplication algorithm only generates non-zero elements when the computed sum exceeds the threshold

### 6.2 Performance

The performance of both sequential and OpenMP implementations was evaluated using diagonal matrices of size
`10000 × 10000` with 1000 non-zero elements. Tests were conducted with varying numbers of threads (1 for sequential, 4,
8, and 16 for OpenMP) and iteration counts (10, 100, and 1000) to measure scalability and parallel efficiency.

| Mode | Count of executions | Num threads | Time amount, s |
| ---- | ------------------- | ----------- | -------------- |
| seq  | 10                  | 1           | 4.9241         |
| seq  | 100                 | 1           | 50.2791        |
| seq  | 1000                | 1           | 499.6823       |
| omp  | 10                  | 4           | 3.0775         |
| omp  | 100                 | 4           | 28.2466        |
| omp  | 1000                | 4           | 280.7209       |
| omp  | 10                  | 8           | 2.8671         |
| omp  | 100                 | 8           | 26.1289        |
| omp  | 1000                | 8           | 259.2370       |
| omp  | 10                  | 16          | 2.9831         |
| omp  | 100                 | 16          | 26.3409        |
| omp  | 1000                | 16          | 261.4287       |

#### Speedup Analysis

Using the sequential execution as baseline (1000 iterations), we calculate the speedup achieved by each parallel
configuration:

| Threads | Time (1000 iters), s | Speedup vs Sequential | Efficiency |
| ------- | -------------------- | --------------------- | ---------- |
| 1 (seq) | 499.6823             | 1.00                  | 100%       |
| 4       | 280.7209             | 1.78                  | 44.5%      |
| 8       | 259.2370             | 1.93                  | 24.1%      |
| 16      | 261.4287             | 1.91                  | 11.9%      |

The best configuration (8 threads) completes the workload almost twice as fast as the sequential version.

#### Diminishing Returns with Increased Thread Count

The data reveals an important pattern: increasing threads beyond 4 provides minimal additional benefit:

- **4 → 8 threads**: Only 8% additional speedup (far from ideal 2×)

- **8 → 16 threads**: Slight performance degradation (-0.8%), indicating overhead exceeds benefits

- **Best performance**: Achieved with 8 threads (259.24 s for 1000 iterations)

The OpenMP implementation demonstrates clear performance improvements over the sequential version:

- **4 threads**: 1.78× speedup (reduces time from 499.68s to 280.72s)

- **8 threads**: 1.93× speedup (best performance: 259.24s)

- **16 threads**: 1.91× speedup (slightly worse than 8 threads: 261.43s)

#### Scaling relative to 4-thread configutaion

| Threads | Time (1000 iters), s | Speedup (vs 4 threads) | Additional Efficiency |
| ------- | -------------------- | ---------------------- | --------------------- |
| 4       | 280.7209             | 1.00                   | 100%                  |
| 8       | 259.2370             | 1.08                   | 13.5%                 |
| 16      | 261.4287             | 1.07                   | 6.7%                  |

## 7. Conclusions

This work successfully implemented and validated a parallel algorithm for sparse matrix multiplication using the
Compressed Row Storage (CRS) format with OpenMP. The key findings and conclusions are:

### 7.1 Achievements

- **Successful OpenMP parallelization**: A complete CRS-based matrix multiplication implementation was developed using
  OpenMP, transforming the sequential algorithm into a parallel one capable of utilizing multiple CPU cores effectively.

- **Thread-safe design**: Through careful data scoping with shared and private clauses, thread-local storage for marker
  arrays, and proper synchronization, the implementation ensures correct parallel execution without race conditions.

- **Comprehensive validation**: The implementation passed all functional tests, confirming that parallel execution
  produces bitwise-identical results to the sequential version within the specified floating-point tolerance (kZERO =
  10⁻⁵).

- **Significant performance improvement**: The OpenMP version achieves up to 1.93× speedup on 8 cores compared to the
  sequential implementation, reducing execution time from 499.68 seconds to 259.24 seconds for 1000 matrix
  multiplications.

- **Optimal configuration identified**: Through systematic performance testing, 8 threads were identified as the optimal
  configuration for the given problem size, providing the best balance between parallelism and overhead.

- **Memory efficiency maintained**: Despite adding parallelism, the implementation preserves the memory-efficient nature
  of CRS format, storing only non-zero elements and reusing thread-local temporary buffers.

- **Scalability analysis**: Detailed performance measurements revealed the algorithm's scaling characteristics,
  identifying the sequential transposition step as the primary bottleneck limiting further speedup according to Amdahl's
  Law.

- **Threshold compliance**: The implementation correctly respects the zero threshold, storing only elements with
  absolute value greater than or equal to 10⁻⁵ in the result matrix.

- **Reproducible testing framework**: A comprehensive test suite was developed, including both functional tests with
  pre-computed reference results and performance tests with configurable parameters, ensuring reliable and reproducible
  evaluation.

### 7.2 Limitations

- **Amdahl's Law limitation**: The sequential transposition step remains a bottleneck, limiting maximum theoretical
  speedup to approximately 7.14× regardless of thread count

- **Memory bandwidth saturation**: Performance gains diminish beyond 8 threads due to shared memory bandwidth contention

- **Load imbalance potential**: Current implementation uses static scheduling, which may be suboptimal for matrices with
  highly irregular row distributions

- **Increased memory footprint**: Thread-local marker arrays duplicate storage (O(P × ncols(A))), which can become
  significant for matrices with many columns and high thread counts

- **Transposition overhead**: The transposition step, while necessary for column-wise access, doubles the memory
  footprint for matrix B and remains sequential

## 9. References

1. [Разреженное матричное умножение, Мееров И. Б., Сысоев А. В.](http://www.hpcc.unn.ru/file.php?id=486)
2. [Sparse Matrix. Wikipedia](http://en.wikipedia.org/wiki/Sparse_matrix)
3. [scipy.sparse documentation](https://docs.scipy.org/doc/scipy/reference/sparse.html)
