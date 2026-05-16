# Linear Image Contrast Stretching — ALL

- Student: Zaharov Gleb Mihajlovič, group 3823Б1ПР4
- Technology: ALL / MPI + OpenMP + `std::thread` + oneTBB
- Variant: 28

## 1. Context

The hybrid version demonstrates hierarchical parallelism: work is distributed
between MPI ranks, while different threading technologies are used inside each
process. This version is needed not only for speedup, but also to demonstrate
the ability to combine several execution models in one task.

## 2. Problem Statement

The task is the same as in SEQ: find global `min/max` over the intensity array
and perform linear contrast stretching. Unlike pure threaded backends, this
implementation must coordinate the result between MPI processes.

## 3. Baseline Algorithm

1. Each MPI rank obtains its index range through `GetRankRange`.
2. Inside its range, a rank searches for local `min/max` using `std::thread`.
3. All ranks merge local values through `MPI_Allreduce` and obtain global `min/max`.
4. The array transformation inside a process is performed by three parts: OpenMP,
`std::thread`, and oneTBB.
5. `MPI_Barrier` is called at the end to synchronize completion of all processes.

## 4. Interprocess Scheme

The range for a rank is computed as follows:

```cpp
return {size * rank_id / ranks_count, size * (rank_id + 1) / ranks_count};
```

This covers the whole input array without overlaps between ranks.

The local minimum and maximum are merged as follows:

```cpp
MPI_Allreduce(&local.min, &global.min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
MPI_Allreduce(&local.max, &global.max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
```

`MPI_Allreduce` is used because global `min/max` values are required by every
rank for identical pixel transformation. If the result were needed only by one
process, `MPI_Reduce` would be sufficient, but this implementation requires the
same scale on all processes.

`MPI_Barrier` at the end of `RunImpl` guarantees that all processes have reached
the end of the computation before the task finishes.

## 5. Intraprocess Scheme

Inside a process, three technologies are used:

- `std::thread` for local `min/max` search in the rank range;
- OpenMP for processing the first third of the array;
- `std::thread` for processing the second third of the array;
- oneTBB for processing the last third of the array.

Key processing split:

```cpp
const std::size_t first_border = input.size() / 3;
const std::size_t second_border = (2 * input.size()) / 3;

StretchOmpRange(input, output, 0, first_border, minmax.min, denom);
StretchStlRange(input, output, first_border, second_border, minmax.min, denom);
StretchTbbRange(input, output, second_border, input.size(), minmax.min, denom);
```

Each part writes to a non-overlapping range of the output array. Therefore,
protection of individual `output` elements is not required between these parts.

## 6. Implementation Details

Implementation files:

- `all/include/ops_all.hpp`
- `all/src/ops_all.cpp`

`GetThreadCount` limits the number of threads by both the range size and
`ppc::util::GetNumThreads()`. This prevents creating unnecessary threads on
small ranges.

`FindLocalMinMaxStl` returns neutral values `{255, 0}` for an empty range. This
is important when the number of ranks is greater than the number of elements:
such a rank must not distort the global reduction.

If the image is constant, parallel copying through TBB is used. If `max > min`,
the transformation is performed using the common formula.

## 7. Correctness Check

Correctness of ALL is checked by comparing the result with SEQ and the reference
function `ReferenceLinContrStr`. Tests should be run with different configurations:

| Configuration                                       | What is checked                                             |
|-----------------------------------------------------|-------------------------------------------------------------|
| `1 × 1`                                             | equivalence to sequential logic                             |
| `1 × N`                                             | intraprocess parallelism without distribution between ranks |
| `2 × N`                                             | correctness of `MPI_Allreduce` and rank-range splitting     |
| `P × N`, where `P` is larger than some small inputs | robustness to empty rank ranges                             |

## 8. Experimental Environment

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
export PPC_NUM_PROC=2
export PPC_NUM_THREADS=4
scripts/run_tests.py --running-type=processes --counts 2 4 8
scripts/run_tests.py --running-type=performance
```

## 9. Results

| Ranks | Threads per rank | Total workers | Mode     | Size      | Time, ms | Speedup vs SEQ | Efficiency |
|------:|-----------------:|--------------:|---------:|-----------|---------:|---------------:|-----------:|
| 1     | 1                | 1             | task     | `1 << 20` | 6.853    | 0.64           | 0.64       |
| 1     | 2                | 2             | task     | `1 << 20` | 3.931    | 1.12           | 0.56       |
| 1     | 4                | 4             | task     | `1 << 20` | 3.283    | 1.34           | 0.33       |
| 1     | 8                | 8             | task     | `1 << 20` | 3.834    | 1.15           | 0.14       |
| 1     | 1                | 1             | pipeline | `1 << 20` | 6.582    | 0.69           | 0.69       |
| 1     | 2                | 2             | pipeline | `1 << 20` | 5.363    | 0.85           | 0.42       |
| 1     | 4                | 4             | pipeline | `1 << 20` | 4.093    | 1.11           | 0.28       |
| 1     | 8                | 8             | pipeline | `1 << 20` | 3.199    | 1.42           | 0.18       |

## 10. Conclusions

In the provided logs, `all` was run as a hybrid backend with one MPI rank: the
separate `*_mpi_*` filter did not match any tests. In this configuration, the
version demonstrates correctness and some speedup as the number of threads
increases, but it remains noticeably slower than OMP/STL/TBB because of
additional overhead from several runtimes and MPI synchronization. A full
evaluation of the interprocess part requires separate runs through the
process/MPI runner with `PPC_NUM_PROC > 1`.
