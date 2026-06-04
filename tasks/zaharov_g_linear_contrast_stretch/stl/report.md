# Linear Image Contrast Stretching — STL

- Student: Zaharov Gleb Mihajlovič, group 3823Б1ПР4
- Technology: STL / `std::thread`
- Variant: 28

## 1. Context

The STL version implements parallelism manually through `std::thread`. Unlike
OpenMP and TBB, this backend explicitly defines work ranges, creates threads,
and waits for their completion.

## 2. Problem Statement

The goal is to obtain the same result as in SEQ: perform linear stretching of
the input `uint8_t` intensities to the range `[0, 255]`. SEQ is used as the
reference.

## 3. Baseline Algorithm

The algorithm consists of two parallel stages:

1. search for local minima and maxima in thread ranges;
2. element-wise transformation of the input into the output.

After the first stage, local `{min, max}` pairs are merged into the global pair.
After the second stage, all threads are finished and the output array is ready.

## 4. Parallelization Scheme

The number of threads is chosen as the minimum of the input size and
`ppc::util::GetNumThreads()`. This prevents creating more threads than there are
elements in the input array.

The range for `thread_index` is computed as follows:

```cpp
begin = input.size() * thread_index / thread_count;
end = input.size() * (thread_index + 1) / thread_count;
```

This splitting covers the whole array without gaps or overlaps.

Key thread-launching fragment:

```cpp
for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
  const std::size_t begin = input.size() * thread_index / thread_count;
  const std::size_t end = input.size() * (thread_index + 1) / thread_count;
  threads.emplace_back([&input, &output, begin, end, minmax]() {
    StretchRange(input, output, begin, end, minmax.min, minmax.max);
  });
}

for (auto &thread : threads) {
  thread.join();
}
```

It is important that `join()` is called only after all threads have been
created. Therefore, work can actually run in parallel instead of being
serialized inside the creation loop.

## 5. Implementation Details

Implementation files:

- `stl/include/ops_stl.hpp`
- `stl/src/ops_stl.cpp`

For `min/max` search, each thread writes its result only to its own
`local_minmax[thread_index]` cell. This eliminates write races. During the
transformation, each thread writes only to its own `output[begin:end)` range.
Therefore, neither `mutex` nor `atomic` is required.

If the image is constant, `StretchRange` copies the corresponding input range to
the output. The copy is also safe because thread ranges do not overlap.

## 6. Correctness Check

The STL version is checked with the same functional cases as the other backends.
Special attention is paid to cases where the array size is not divisible by the
number of threads: the integer-division splitting formulas must cover the tail
of the array in the last ranges.

The result must match `ReferenceLinContrStr` for `PPC_NUM_THREADS = 1`, `2`,
`4`, and other values allowed by the course infrastructure.

## 7. Experimental Environment

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
export PPC_NUM_THREADS=4
scripts/run_tests.py --running-type=threads --counts 1 2 4 8
scripts/run_tests.py --running-type=performance
```

## 8. Results

| Threads | Mode     | Size      | Time, ms | Speedup vs SEQ | Efficiency |
|--------:|---------:|-----------|---------:|---------------:|-----------:|
| 1       | task     | `1 << 20` | 4.674    | 0.94           | 0.94       |
| 2       | task     | `1 << 20` | 1.701    | 2.58           | 1.29       |
| 4       | task     | `1 << 20` | 1.051    | 4.18           | 1.04       |
| 8       | task     | `1 << 20` | 1.054    | 4.17           | 0.52       |
| 1       | pipeline | `1 << 20` | 4.454    | 1.02           | 1.02       |
| 2       | pipeline | `1 << 20` | 1.805    | 2.51           | 1.26       |
| 4       | pipeline | `1 << 20` | 1.039    | 4.37           | 1.09       |
| 8       | pipeline | `1 << 20` | 1.108    | 4.10           | 0.51       |

## 9. Conclusions

According to the measurements, the STL version is most efficient with 4 threads:
`1.051 ms` for `task` and `1.039 ms` for `pipeline`. With 8 threads there is no
additional gain, because thread creation/waiting costs and memory pressure start
to dominate.
