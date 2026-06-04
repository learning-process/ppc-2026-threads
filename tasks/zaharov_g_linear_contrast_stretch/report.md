# Linear Image Contrast Stretching

- Student: Zaharov Gleb Mihajlovič, group 3823Б1ПР4
- Variant: 28
- Task directory: `tasks/zaharov_g_linear_contrast_stretch`
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`,
`stl/report.md`, `all/report.md`

## 1. Introduction

The task is to perform linear contrast stretching for a grayscale image. The
image is represented as a one-dimensional byte array `uint8_t`, where each
element is a pixel intensity from `0` to `255`. The algorithm first finds the
minimum and maximum intensity, then maps every pixel to the new range `[0, 255]`
using a linear formula. This task is well suited for comparing parallel
programming models: finding the minimum and maximum is a reduction, while
transforming pixels is an independent element-wise pass over the array.

The project implements five backends: the sequential `seq` version, the OpenMP
`omp` version, the oneTBB `tbb` version, the manual thread-based `stl` version,
and the hybrid `all` version combining MPI, OpenMP, `std::thread`, and oneTBB.

## 2. Unified Problem Statement

Input data: a non-empty array `InType = std::vector<uint8_t>`.

Output data: an array `OutType = std::vector<uint8_t>` of the same length,
containing pixels after linear contrast stretching.

For the input array `input`, the following values are defined:

```text
min_el = min(input)
max_el = max(input)
```

If `max_el > min_el`, every pixel is transformed as follows:

```text
output[i] = (input[i] - min_el) * 255 / (max_el - min_el)
```

If all pixels are equal, then `max_el == min_el`; in this case contrast
stretching is impossible, so the result is a copy of the input image. This
prevents division by zero and preserves correct behavior for constant images.

Correctness criterion: the result of every parallel implementation must be
byte-for-byte equal to the reference sequential formula for all functional and
performance tests.

## 3. Unified Experimental Methodology

All implementations use the common data types from `common/include/common.hpp`
and the same task lifecycle:

1. `ValidationImpl` checks that the input array is non-empty.
2. `PreProcessingImpl` allocates an output array of the required size.
3. `RunImpl` performs the `min/max` search and pixel transformation.
4. `PostProcessingImpl` checks that the output array is not empty.

Functional tests are located in `tests/functional/main.cpp` and run all enabled
backends through the common `BaseRunFuncTests` mechanism. Four input types are
checked: a short shifted range, a constant image, the full `[0, 255]` range, and
pseudo-random data.

Performance tests are located in `tests/performance/main.cpp`. They use an array
of size `1 << 20` with a fixed generator seed, so the input is reproducible. All
backends are compared against the same reference result.

Performance analysis formulas:

```text
speedup = T_seq / T_parallel
efficiency = speedup / workers
```

For `seq`, the number of workers is `1`. For `omp`, `tbb`, and `stl`, workers
are threads. For `all`, the configuration must be written as `ranks × threads`,
and efficiency is normalized by `total_workers = ranks * threads_per_rank`.

## 4. Experimental Environment

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

Basic build and run commands:

```bash
git submodule update --init --recursive --depth=1

cmake -S . -B build \
  -D USE_FUNC_TESTS=ON \
  -D USE_PERF_TESTS=ON \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel

export PPC_NUM_THREADS=4
scripts/run_tests.py --running-type=threads --counts 1 2 4 8

export PPC_NUM_PROC=2
scripts/run_tests.py --running-type=processes --counts 2 4

scripts/run_tests.py --running-type=performance
```

Nix shell I used:

```nix
{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  nativeBuildInputs = with pkgs; [
    gcc
    gdb
    cmake
    pkg-config
    ninja
    openmpi
    clang-tools
    lldb_21
  ];

  buildInputs = with pkgs; [
    openmpi
    valgrind
    python3Packages.python
    python3Packages.virtualenv
  ];

  shellHook = ''
    export CC=${pkgs.gcc}/bin/gcc
    export CXX=${pkgs.gcc}/bin/g++
    export OMPI_CC=${pkgs.gcc}/bin/gcc
    export OMPI_CXX=${pkgs.gcc}/bin/g++

    # Пути к локальным библиотекам (после сборки)
    export REPO="/mnt/D/Coding/UNN/ppc-2026-threads"
    export LD_LIBRARY_PATH="$REPO/build/ppc_googletest/install/lib64:$REPO/build/ppc_libenvpp/install/lib64:$REPO/build/ppc_onetbb/install/lib64:$LD_LIBRARY_PATH"
    export LIBRARY_PATH="$LD_LIBRARY_PATH"
    export CPLUS_INCLUDE_PATH="$REPO/build/ppc_googletest/install/include:$REPO/build/ppc_libenvpp/install/include:$REPO/build/ppc_onetbb/install/include:$CPLUS_INCLUDE_PATH"

    if [ ! -d ".venv" ]; then
      python -m venv .venv
    fi

    source .venv/bin/activate
  '';
}
```

When measurements are performed on Linux, background processes should be
minimized and the same thread configuration should be used for all runs. If the
CPU governor or CPU affinity was not fixed, this must be stated explicitly when
interpreting the results.

## 5. Correctness Summary

All implementations are compared with the same reference function
`ReferenceLinContrStr` defined in the tests. It repeats the mathematical
statement of the task and does not depend on any particular backend.

Tested cases:

| Case         | Purpose                               | Expected property                                                |
|--------------|---------------------------------------|------------------------------------------------------------------|
| `shifted`    | small intensity range `{50, 75, 100}` | the result is stretched to `{0, 127, 255}` with integer division |
| `constant`   | all pixels are equal                  | output equals input                                              |
| `full_range` | values from `0` to `255`              | output equals input                                              |
| `random`     | medium pseudo-random input            | result matches the reference                                     |

There is also a separate check that `ValidationImpl` rejects an empty input.

## 6. Aggregated Results

The final results below are based on the provided performance logs. Time is
reported in milliseconds. For SEQ, the median of all provided SEQ runs is used:
`task = 4.392 ms`, `pipeline = 4.539 ms`. For the other backends, rows
correspond to run groups interpreted as `PPC_NUM_THREADS = 1, 2, 4, 8`.

### Performance mode: task

| Backend | Workers / configuration | Median time, ms | Speedup vs SEQ | Efficiency | Notes                                             |
|---------|------------------------:|----------------:|---------------:|-----------:|---------------------------------------------------|
| SEQ     |                   1     |           4.392 |           1.00 |       1.00 | median over all SEQ task runs                     |
| OMP     |                   1     |           2.888 |           1.52 |       1.52 | `min/max` + stretch in OpenMP                     |
| OMP     |                   2     |           1.567 |           2.80 |       1.40 | best balance without strong oversubscription      |
| OMP     |                   4     |           0.934 |           4.70 |       1.18 | best OMP task result                              |
| OMP     |                   8     |           0.944 |           4.65 |       0.58 | speedup almost stops growing after 4 threads      |
| STL     |                   1     |           4.674 |           0.94 |       0.94 | manual thread-based version is close to SEQ       |
| STL     |                   2     |           1.701 |           2.58 |       1.29 | significant gain with two threads                 |
| STL     |                   4     |           1.051 |           4.18 |       1.04 | best STL task result                              |
| STL     |                   8     |           1.054 |           4.17 |       0.52 | saturation after 4 threads                        |
| TBB     |                   1     |           2.638 |           1.66 |       1.66 | runtime provides a fast pass even with one worker |
| TBB     |                   2     |           1.386 |           3.17 |       1.58 | stable gain                                       |
| TBB     |                   4     |           0.740 |           5.93 |       1.48 | best task result among pure threaded backends     |
| TBB     |                   8     |           0.910 |           4.82 |       0.60 | worse than with 4 workers                         |
| ALL     |                   1 × 1 |           6.853 |           0.64 |       0.64 | hybrid overhead is higher than SEQ                |
| ALL     |                   1 × 2 |           3.931 |           1.12 |       0.56 | speedup is limited by overhead                    |
| ALL     |                   1 × 4 |           3.283 |           1.34 |       0.33 | best ALL task result                              |
| ALL     |                   1 × 8 |           3.834 |           1.15 |       0.14 | runtime oversubscription worsens the result       |

### Performance mode: pipeline

| Backend | Workers / configuration | Median time, ms | Speedup vs SEQ | Efficiency | Notes                                                  |
|---------|------------------------:|----------------:|---------------:|-----------:|--------------------------------------------------------|
| SEQ     | 1                       | 4.539           | 1.00           | 1.00       | median over all SEQ pipeline runs                      |
| OMP     | 1                       | 2.904           | 1.56           | 1.56       | `min/max` + stretch in OpenMP                          |
| OMP     | 2                       | 1.811           | 2.51           | 1.25       | stable speedup                                         |
| OMP     | 4                       | 0.966           | 4.70           | 1.17       | best OMP pipeline result                               |
| OMP     | 8                       | 1.031           | 4.40           | 0.55       | speedup does not improve after 4 threads               |
| STL     | 1                       | 4.454           | 1.02           | 1.02       | close to SEQ                                           |
| STL     | 2                       | 1.805           | 2.51           | 1.26       | significant gain                                       |
| STL     | 4                       | 1.039           | 4.37           | 1.09       | best STL pipeline result                               |
| STL     | 8                       | 1.108           | 4.10           | 0.51       | saturation after 4 threads                             |
| TBB     | 1                       | 2.966           | 1.53           | 1.53       | runtime overhead is compensated by efficient splitting |
| TBB     | 2                       | 1.652           | 2.75           | 1.37       | stable gain                                            |
| TBB     | 4                       | 1.024           | 4.43           | 1.11       | close to OMP/STL                                       |
| TBB     | 8                       | 0.727           | 6.24           | 0.78       | best pipeline result                                   |
| ALL     | 1 × 1                   | 6.582           | 0.69           | 0.69       | hybrid overhead is higher than SEQ                     |
| ALL     | 1 × 2                   | 5.363           | 0.85           | 0.42       | still slower than SEQ                                  |
| ALL     | 1 × 4                   | 4.093           | 1.11           | 0.28       | limited speedup                                        |
| ALL     | 1 × 8                   | 3.199           | 1.42           | 0.18       | best ALL pipeline result                               |

## 7. Interpretation of Differences

- `SEQ` is the correctness reference and the baseline for computing speedup. Its
computational complexity is linear in the number of pixels, and memory
consumption is also linear because of the output array.

- `OMP` parallelizes both computational phases. The first phase finds local
minima and maxima in threads and then combines them in one thread. The second
phase transforms independent elements of the output array in parallel. The main
overhead comes from creating parallel regions and merging local results.

- `STL` uses explicit splitting of the array into contiguous ranges. Each thread
writes only to its own part of the output array, so the main phase does not
require a `mutex` or `atomic`. The main overhead is manual thread creation and
waiting with `join`.

- `TBB` delegates range splitting and task scheduling to the oneTBB runtime.
`parallel_reduce` is used for finding the minimum and maximum, and
`parallel_for` is used for the transformation. This version is convenient
because it does not require manual thread management, but runtime overhead may
be noticeable for small inputs.

- `ALL` uses a hierarchical scheme. First, MPI ranks receive their input ranges
and find global `min/max` values with `MPI_Allreduce`. Then transformation
inside a process is split between OpenMP, `std::thread`, and oneTBB. This
version demonstrates the hybrid model, but its speedup is limited by MPI
synchronization and the additional overhead of several runtimes.

## 8. Conclusion

For linear contrast stretching, the expected main benefit comes from backends
that efficiently parallelize two linear passes over the array. The task has
almost no dependencies between elements, so parallel versions can work correctly
without synchronization on every pixel. The most important performance factors
are the input size, the number of threads, the cost of the `min/max` reduction,
runtime overhead, and memory bandwidth.

According to the provided measurements, the best version in `task` mode is TBB
with 4 workers (`0.740 ms`), while the best version in `pipeline` mode is TBB
with 8 workers (`0.727 ms`). OMP and STL also provide stable speedup up to 4
threads, after which the improvement saturates. The hybrid `all` version with
one MPI rank is slower than the pure threaded backends because of additional
runtime and synchronization overhead.

## 9. References

- Course documentation and the `example_threads` structure.
- OpenMP specification.
- oneTBB documentation.
- MPI Forum documentation.
- C++ reference documentation for `std::thread`.
