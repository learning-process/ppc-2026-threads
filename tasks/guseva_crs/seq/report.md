# Sparse Matrix Multiplication. Double Precision Elements. Compressed Row Storage (CRS) Format

- _Student_: Гусева Алёна Сергеевна, group 3823Б1ФИ2
- _Technology_: SEQ
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

### 1.3 Expected Outcome

The expected outcome of this work is a sequential implementation of sparse matrix multiplication for matrices stored in
CRS format. The implementation will read matrices from test files, perform multiplication with correctness verification,
and produce the result matrix also in CRS format.

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

The baseline algorithm for sparse matrix multiplication in CRS format consists of several steps:

### 3.1 Matrix Transposition

Before performing the multiplication, matrix B is transposed. Transposition converts the column-oriented access pattern
into a row-oriented one, which significantly improves cache efficiency and simplifies the implementation. The
transposition algorithm works as follows:

1. Create a new CRS structure for the transposed matrix Bᵀ

2. Set `nrows = B.ncols`, `ncols = B.nrows`

3. Count the number of elements in each row of `Bᵀ` by iterating through `B.cols` and incrementing counters

4. Build the row pointers for `Bᵀ` using prefix sums

5. Fill the values and column indices of `Bᵀ` by iterating through the original matrix `B` and placing each element in
   its new position

### 3.2 Multiplication Algorithm

Once B is transposed, the multiplication proceeds row by row: For each row i of matrix A:

1. Initialize a temporary array or marker to track which columns of `A` have non-zero elements in this row

2. For each non-zero element in row `i` of `A`, record its column index and position

3. For each row `j` of `Bᵀ` (which corresponds to column `j` of `B`):
    - Initialize `sum = 0`

    - For each non-zero element in row `j` of `Bᵀ`:
        - If the column index of this element (which corresponds to a row index in original `B`) matches a recorded
          column from row `i` of `A`

        - Multiply the corresponding values from `A` and `Bᵀ` and add to `sum`

    - If `|sum| > 10⁻⁵`, add this element to the result matrix `C` at position `(i, j)`

### 3.3 Complexity Analysis

**Time complexity:** `O(nz(A) + nz(B) + number_of_non-zero_products)` for the actual multiplication, plus `O(nz(B))` for
transposition

**Space complexity**: `O(nz(A) + nz(B) + nz(C))` for storing matrices, plus temporary arrays of size `O(ncols(A))` for
each row processing

## 4. Implementation Details

### 4.1 Code Structure

```└── tasks/
    └── guseva_crs/
        ├── common/
        |   ├── common.hpp  . . . . . . . . . . . . Common type definitions
        |   |                                       and constants, CRS struct definition
        |   ├── multiplier.hpp  . . . . . . . . . . Base class Multiplier
        |   └── test_reader.hpp . . . . . . . . . . Functions to read test data and store them in CRS format
        |
        ├── data/ . . . . . . . . . . . . . . . . . Directory for test files in .txt format
        |
        ├── seq/
        |   ├── inlude/
        |   |   ├── ops_seq.hpp . . . . . . . . . . Sequential task class declaration
        |   |   └── multiplier_seq.hpp  . . . . . . Realisation on sequentional multiplier
                                                    inherited from Multiplier defined in ./common/multiplier.hpp
        |   |
        |   └── src/
        |       └── ops_seq.cpp . . . . . . . . . . Sequential implementation
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
- `class Multiplier` inherited from `class Multiplier` to perform dot product sequentionally. Contains all the
  executions
- `class GusevaCRSMatMulSeq` implements the repo's `BaseTask` interface for correct task verification and acceptance
- `class GusevaMatMulCRSFuncTest` to produce functional tests
- `class GusevaMatMulCRSPerfTest` to produce performance tests

### 4.3 Important Assumptions And Corner Cases

**Assumptions:**

- Matrices are valid and properly formatted in CRS with 0-based indexing

- Row pointers satisfy: row_ptrs[0] = 0, row_ptrs[nrows] = nz

- Column indices for each row are sorted (typical CRS requirement)

- Input matrices contain only double-precision floating-point values

**Corner Cases Handled:**

- Empty matrices: Matrices with nz = 0 are handled correctly

- Zero threshold: Values below 10⁻⁵ are considered zero and not stored in result

- Invalid dimensions: Validation checks that input matrices are multiplicatable

- Single-element matrices: Multiplication works for 1×1 matrices

- Rectangular matrices: Handles non-square matrices correctly

### 4.4 Memory Usage Considerations

The implementation is designed with memory efficiency in mind:

1. Sparse storage: Only non-zero elements are stored, saving significant memory for sparse matrices

2. In-place transposition: The transpose operation creates a new matrix but avoids unnecessary copying

3. Temporary arrays: The temp vector is reused for each row, sized to the number of columns in A

4. Result construction: The result matrix is built incrementally without pre-allocation of full dense storage

**Memory complexity:**

1. Input matrices: O(nz(A) + nz(B)) + O(nrows(A) + nrows(B) + ncols(A) + ncols(B)) for index arrays

2. Transposed matrix: O(nz(B)) + O(ncols(B)) for row pointers

3. Temporary storage: O(ncols(A)) for the marker array

4. Output matrix: O(nz(C)) + O(nrows(C)) for row pointers

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

- **Environment**: Not needed for sequentional setup
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

The correctness of the sequential sparse matrix multiplication implementation was verified through a comprehensive
testing approach involving multiple verification strategies:

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

The `Equal` function implements a robust comparison that checks:

- Structural equality: same number of non-zero elements, same dimensions

- Value equality: each corresponding value pair must satisfy `|a - b| < kZERO`

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

The performance of the sequential implementation was evaluated using test data with varying iteration counts. The
measurements show the execution time for multiplying sparse matrices of fixed size across multiple runs:

| Mode | Count | Time, s  |
| ---- | ----- | -------- |
| seq  | 10    | 4.9241   |
| seq  | 100   | 50.2791  |
| seq  | 1000  | 499.6823 |

The performance data reveals approximately linear scaling with the number of multiplications performed. The time
increases from 14.33 seconds for 10 iterations to 1405.68 seconds for 1000 iterations, showing a growth factor of about
98× for a 100× increase in workload. This near-linear scaling is expected for a sequential implementation as there is no
parallelization overhead.

However, the implementation exhibits inherent scalability limitations due to its sequential nature. The execution time
grows proportionally with the number of non-zero elements in the input matrices and the density of the resulting matrix.
For matrices with high sparsity ratios, the algorithm performs efficiently as it only processes non-zero elements.
Conversely, as matrix density increases toward full matrices, performance degrades significantly since the algorithm
would need to handle O(n³) operations despite still using sparse storage.

The primary bottleneck in the current implementation is the row-by-row processing with temporary array initialization
for each row of matrix A. For matrices with many columns, this repeated initialization contributes noticeable overhead.
Additionally, the transposition step, while enabling more efficient column access, doubles the memory footprint for
matrix B and adds an extra pass through all non-zero elements.

## 7. Conclusions

This work successfully implemented and validated a sequential algorithm for sparse matrix multiplication using the
Compressed Row Storage (CRS) format. The key findings and conclusions are:

### 7.1 Achievements

- A complete CRS-based matrix multiplication implementation was developed with proper handling of sparse data structures

- The algorithm correctly handles matrix transposition to optimize column-wise access patterns

- Comprehensive validation ensures correctness through reference comparisons and invariant checking

- The implementation respects the zero threshold (10⁻⁵), storing only elements with significant magnitude

- Memory efficiency is achieved by storing only non-zero elements and reusing temporary buffers

### 7.2 Limitations

- The sequential nature limits performance on large matrices, with execution time scaling linearly with workload

- For dense matrices, the algorithm loses efficiency as the sparse representation becomes less beneficial

- The temporary marker array, sized to the number of columns in A, can become a memory bottleneck for matrices with many
  columns

- Floating-point accumulation may lead to cancellation errors, though the threshold helps mitigate this

- Single-threaded execution cannot utilize modern multi-core processors for performance gains

The implementation provides a solid foundation for sparse matrix multiplication that balances correctness, memory
efficiency, and reasonable performance for moderately sized problems. However, for large-scale applications or real-time
requirements, parallelization strategies would be necessary to overcome the inherent scalability limitations of the
sequential approach.

## 9. References

1. [Разреженное матричное умножение, Мееров И. Б., Сысоев А. В.](http://www.hpcc.unn.ru/file.php?id=486)
2. [Sparse Matrix. Wikipedia](http://en.wikipedia.org/wiki/Sparse_matrix)
3. [scipy.sparse documentation](https://docs.scipy.org/doc/scipy/reference/sparse.html)

## Appendix

### 1. `generate_tests.py`

```python
import numpy as np
from scipy.sparse import csr_matrix as crs, random

def generate_matrix(rows, cols, density):
    return random(rows, cols, density=density, format='csr', data_rvs=lambda s: np.random.randint(1, 100, s) / 2)

def generate_tests(params):
    for size_a, size_b, density, name in params:
        a = generate_matrix(*size_a, density[0])
        b = generate_matrix(*size_b, density[1])
        a.sort_indices()
        b.sort_indices()
        c = a.dot(b)
        c.sort_indices()
        drop_test2file(name + '.txt', a, b, c)

def write_matrix(file, a):
    file.write(f"{a.nnz} {a.shape[0]} {a.shape[1]}\n")
    [file.write(f"{str(x)} ") for x in a.data]
    file.write('\n')
    [file.write(f"{str(x)} ") for x in a.indices]
    file.write('\n')
    [file.write(f"{str(x)} ") for x in a.indptr]
    file.write('\n')


def drop_test2file(filename, a, b, c):
    with open(filename, 'w', encoding='utf-8') as file:
        write_matrix(file, a)
        write_matrix(file, b)
        write_matrix(file, c)



if __name__ == '__main__':
    generate_tests([
        ((5, 5), (5, 5), (0.2, 1), 'sparse_dense'),
        ((10, 10), (10, 10), (1, 0.2), 'dense_sparse'),
        ((15, 15), (15, 15), (0.1, 0.1), 'double_sparse1'),
        ((13, 13), (13, 13), (0.1, 0.1), 'double_sparse2'),
        ((23, 23), (23, 23), (0.1, 0.1), 'double_sparse3'),
        ((31, 31), (31, 31), (0.1, 0.1), 'double_sparse4'),
        ((5e3, 5e3), (5e3, 5e3), (.005, .005), 'perf')
    ]
    )


```
