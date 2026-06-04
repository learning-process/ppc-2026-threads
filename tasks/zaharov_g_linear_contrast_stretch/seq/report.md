# Linear Image Contrast Stretching — SEQ

- Student: Zaharov Gleb Mihajlovič, group 3823Б1ПР4
- Technology: SEQ
- Variant: 28

## 1. Context

The sequential version is the correctness reference for all other backends. It
implements the mathematical statement directly: first it finds the minimum and
maximum intensity, then it applies a linear transformation to every pixel.

## 2. Problem Statement

Input: a non-empty `std::vector<uint8_t>` representing pixel intensities.

Output: an array of the same length. If the input contains at least two
different values, the minimum intensity is mapped to `0`, the maximum intensity
is mapped to `255`, and all other values are scaled linearly. If all elements
are equal, the output is a copy of the input.

An empty input is rejected in `ValidationImpl`, because its minimum and maximum
cannot be defined.

## 3. Baseline Algorithm

1. Get references to the input and output arrays.
2. Find `min_el` and `max_el` with `std::ranges::minmax_element`.
3. If `max_el > min_el`, iterate over all pixels and compute the new value.
4. If `max_el == min_el`, copy the input to the output.

Time complexity: `O(n)`, where `n` is the number of pixels. The `min/max` search
requires one linear pass, and the transformation requires a second linear pass.

Additional memory: `O(n)` for the output array. Auxiliary variables take `O(1)`
memory.

## 4. Implementation Details

Implementation files:

- `seq/include/ops_seq.hpp`
- `seq/src/ops_seq.cpp`

Key fragment:

```cpp
std::ranges::minmax_element(input);

for (size_t i = 0; i < input.size(); ++i) {
  int value = (static_cast<int>(input[i]) - min_el) * 255 / denom;
  output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
}
```

`ValidationImpl` returns `false` for an empty input. `PreProcessingImpl` resizes
the output array to the input size in advance. `RunImpl` performs the main
computation. `PostProcessingImpl` checks that the output array is not empty.

Using `int` in the formula avoids overflow of `uint8_t` during multiplication.
`std::clamp` additionally fixes the result range, although mathematically the
value should already be within `[0, 255]` when `min/max` are correct.

## 5. Correctness Check

SEQ is compared against the reference function `ReferenceLinContrStr` from the
tests. The following cases are checked:

| Input                   | What is checked | Expected result          |
|-------------------------|-----------------|--------------------------|
| `{50, 75, 100}`         | shifted range   | `{0, 127, 255}`          |
| ten values equal to `7` | constant image  | copy of the input        |
| `0..255`                | full range      | copy of the input        |
| pseudo-random array     | general case    | match with the reference |

There is also the `ValidationRejectsEmptyInput` test, which checks that an empty array is rejected.

## 6. Experimental Environment

Final measurements were performed in the following environment:

**CPU**:                    Intel Core i5-10600KF (6 cores, 12 threads, 4.1GHz base)
**RAM**:                    32 GB
**OS**:                     NixOS
**Compiler**:               GCC 15.2.0
**CMake**:                  4.1.2
**Ninja**:                  1.13.2
**CMake build type**:       `Release`
**C++ standard**:           C++23
**`PPC_NUM_THREADS`**:      1, 2, 4, 8
**Performance input size**: `1 << 20` bytes

Commands:

```bash
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel
scripts/run_tests.py --running-type=threads --counts 1
scripts/run_tests.py --running-type=performance
```

## 7. Results

| Mode     | Size      | Workers | Time, ms | Comment                                              |
|----------|----------:|--------:|---------:|------------------------------------------------------|
| task     | `1 << 20` | 1       | 4.392    | baseline; median over all provided SEQ task runs     |
| pipeline | `1 << 20` | 1       | 4.539    | baseline; median over all provided SEQ pipeline runs |

The measured SEQ time is used as the denominator when computing speedup for the
other implementations.

## 8. Conclusions

The SEQ version is simple, deterministic, and contains no synchronization. It
defines the reference result and the performance baseline. All parallel
implementations must match it byte-for-byte in functional and performance tests.
