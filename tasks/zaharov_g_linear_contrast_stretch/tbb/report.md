# Linear Image Contrast Stretching — TBB

- Student: Zaharov Gleb Mihajlovič, group 3823Б1ПР4
- Technology: oneTBB
- Variant: 28

## 1. Context

The TBB version uses the task-oriented oneTBB model. The programmer describes
the work range and the operation over the range, while the runtime schedules
tasks between worker threads.

## 2. Problem Statement

The statement is the same as in SEQ: find the minimum and maximum intensity,
then linearly scale every pixel to the `[0, 255]` range. The sequential version
is the correctness reference.

## 3. Baseline Algorithm

The algorithm is split into two oneTBB operations:

1. `parallel_reduce` to find global `min/max`;
2. `parallel_for` to transform pixel ranges.

## 4. Parallelization Scheme

The reduction uses `oneapi::tbb::parallel_reduce` over
`blocked_range<std::size_t>(0, input.size())`. Each subrange returns a local
`{min, max}` pair, and then pairs are merged by `MergeMinMax`.

Key fragment:

```cpp
return oneapi::tbb::parallel_reduce(
    oneapi::tbb::blocked_range<std::size_t>(0, input.size()),
    MinMax{std::numeric_limits<uint8_t>::max(), std::numeric_limits<uint8_t>::min()},
    [&input](const oneapi::tbb::blocked_range<std::size_t> &range, MinMax current) {
      for (std::size_t i = range.begin(); i != range.end(); ++i) {
        const int value = static_cast<int>(input[i]);
        current.min = std::min(current.min, value);
        current.max = std::max(current.max, value);
      }
      return current;
    },
    MergeMinMax);
```

The transformation uses `parallel_for` over the same range type. Each subrange
writes to an independent part of the output array, so element-level
synchronization is not required.

No explicit `grainsize` is specified, so oneTBB chooses splitting automatically.
This simplifies the code, but actual performance depends on the runtime strategy
and input size.

## 5. Implementation Details

Implementation files:

- `tbb/include/ops_tbb.hpp`
- `tbb/src/ops_tbb.cpp`

In `FindMinMax`, the local result is passed by value, so different tasks do not
write to a shared accumulator. Merging is performed through the pure function
`MergeMinMax`.

`StretchImage` distinguishes two cases. If `max > min`, every pixel is scaled
using the formula. If the image is constant, the input is copied to the output
in parallel.

## 6. Correctness Check

The TBB version is compared with the reference function in all functional tests.
Since the order of subrange processing in TBB can differ between runs, it is
important that the reduction operation is associative for the `min/max` pair and
that the element-wise transformation does not depend on iteration order.

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

| Workers | Mode     | Size      | Time, ms | Speedup vs SEQ | Efficiency |
|--------:|---------:|-----------|---------:|---------------:|-----------:|
| 1       | task     | `1 << 20` | 2.638    | 1.66           | 1.66       |
| 2       | task     | `1 << 20` | 1.386    | 3.17           | 1.58       |
| 4       | task     | `1 << 20` | 0.740    | 5.93           | 1.48       |
| 8       | task     | `1 << 20` | 0.910    | 4.82           | 0.60       |
| 1       | pipeline | `1 << 20` | 2.966    | 1.53           | 1.53       |
| 2       | pipeline | `1 << 20` | 1.652    | 2.75           | 1.37       |
| 4       | pipeline | `1 << 20` | 1.024    | 4.43           | 1.11       |
| 8       | pipeline | `1 << 20` | 0.727    | 6.24           | 0.78       |

## 9. Conclusions

According to the measurements, TBB showed the best `task` result with 4 workers
(`0.740 ms`) and the best `pipeline` result with 8 workers (`0.727 ms`). For
this task, TBB schedules the two linear passes efficiently, but the 8-worker
`task` result is worse than the 4-worker result, which indicates
memory-bandwidth saturation.
