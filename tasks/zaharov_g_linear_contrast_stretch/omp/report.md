# Linear Image Contrast Stretching — OMP

- Student: Zaharov Gleb Mihajlovič, group 3823Б1ПР4
- Technology: OMP
- Variant: 28

## 1. Context

The OpenMP version maps the sequential algorithm to the fork-join model. The
task contains two naturally parallel parts: a reduction for finding the minimum
and maximum, and an independent transformation of every pixel.

## 2. Problem Statement

The statement is the same as in SEQ: the input array of `uint8_t` intensities
must be transformed so that its intensity range is stretched to `[0, 255]`. SEQ
is used as the baseline and the correctness reference.

## 3. Baseline Algorithm

The algorithm keeps two phases:

1. `min/max` search;
2. transformation of all elements using the linear stretching formula.

In the OpenMP implementation, both phases are parallelized.

## 4. Parallelization Scheme

The `FindMinMax` function creates a parallel region with `ppc::util::GetNumThreads()` threads. Each thread
stores a local `{min, max}` pair and processes its part of the input array.

Key fragment:

```cpp
#pragma omp parallel default(none) shared(input, local_minmax, size) num_threads(thread_count)
{
  MinMax current{std::numeric_limits<uint8_t>::max(), std::numeric_limits<uint8_t>::min()};
#pragma omp for nowait
  for (std::int64_t i = 0; i < size; ++i) {
    const int value = static_cast<int>(input[static_cast<std::size_t>(i)]);
    current.min = std::min(current.min, value);
    current.max = std::max(current.max, value);
  }
  local_minmax[static_cast<std::size_t>(omp_get_thread_num())] = current;
}
```

`input`, `local_minmax`, and `size` are `shared`. The `current` variable is
declared inside the parallel region, so it is local to each thread. The
`default(none)` directive forces explicit data-sharing attributes and reduces
the risk of hidden mistakes.

After the parallel region finishes, local pairs are merged sequentially. This is
safe because all threads have already finished writing to `local_minmax`.

The second phase uses `#pragma omp parallel for`: each thread writes only to its
own output indices. There are no data races, because each `output[i]` element is
computed by exactly one iteration.

## 5. Implementation Details

Implementation files:

- `omp/include/ops_omp.hpp`
- `omp/src/ops_omp.cpp`

When `max > min`, parallel scaling is performed. For a constant image, the input
is copied to the output in parallel.

No explicit `schedule` is specified, so the default behavior of the OpenMP
runtime is used. For this task, the iteration workload is uniform, so static
distribution or the runtime default policy should produce similar behavior.

## 6. Correctness Check

Correctness is checked by comparing with SEQ/`ReferenceLinContrStr` for all
functional cases: `shifted`, `constant`, `full_range`, and `random`. Tests
should also be run with different `PPC_NUM_THREADS` values, for example `1`,
`2`, and `4`, to ensure that the result does not depend on the number of threads.

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
| 1       | task     | `1 << 20` | 2.888    | 1.52           | 1.52       |
| 2       | task     | `1 << 20` | 1.567    | 2.80           | 1.40       |
| 4       | task     | `1 << 20` | 0.934    | 4.70           | 1.18       |
| 8       | task     | `1 << 20` | 0.944    | 4.65           | 0.58       |
| 1       | pipeline | `1 << 20` | 2.904    | 1.56           | 1.56       |
| 2       | pipeline | `1 << 20` | 1.811    | 2.51           | 1.25       |
| 4       | pipeline | `1 << 20` | 0.966    | 4.70           | 1.17       |
| 8       | pipeline | `1 << 20` | 1.031    | 4.40           | 0.55       |

## 9. Conclusions

According to the measurements, OMP performs best with 4 threads: `0.934 ms` in
`task` mode and `0.966 ms` in `pipeline` mode. Increasing the thread count to 8
almost does not improve speedup, which indicates memory-bandwidth saturation and
overhead from parallel regions.
