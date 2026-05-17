# Сводный отчет: сортировка Хоара с простым слиянием

## Введение и постановка

Работа сравнивает реализации сортировки `std::vector<int>` по неубыванию для backend-ов `seq`, `omp`, `tbb`, `stl` и
`all`. Единые типы входа/выхода заданы в [`common/include/common.hpp`](common/include/common.hpp#L11). Корректность
определяется совпадением с результатом `std::ranges::sort` в функциональном тесте
([`tests/functional/main.cpp`](tests/functional/main.cpp#L35)).

## Единая методика эксперимента

Performance-тест генерирует `N=100000` случайных `int` из диапазона `[-1000000, 1000000]`
([`tests/performance/main.cpp`](tests/performance/main.cpp#L20)). Использован `TaskRun`; режим выбирается в раннере
([`modules/util/include/perf_test_util.hpp`](../../modules/util/include/perf_test_util.hpp#L87)). Speedup считался как
`T_seq / T_backend`, efficiency как `speedup / workers`. Baseline: `seq`, `PPC_NUM_THREADS=1`,
`T_seq = 0.0058254364 s`. Для OMP/TBB применялись `PPC_NUM_THREADS=1,2,4`; TBB дополнительно ограничивается
`tbb::global_control` ([`modules/runners/src/runners.cpp`](../../modules/runners/src/runners.cpp#L150)). STL и локальная
часть ALL берут число потоков из `std::thread::hardware_concurrency()`; на этой машине `nproc = 12`, поэтому такие
строки обозначены как auto workers.
Число повторов по умолчанию равно 5 ([`modules/performance/include/performance.hpp`](../../modules/performance/include/performance.hpp#L21)).
Ограничение методики: `TaskRun` повторяет только `Run()` после одного `PreProcessing()`
([`modules/performance/include/performance.hpp`](../../modules/performance/include/performance.hpp#L62)), а вход
создается заново для каждого gtest-параметра через `std::random_device` ([`tests/performance/main.cpp`](tests/performance/main.cpp#L24)).
Поэтому дополнительно сняты pipeline-замеры, где каждый повтор проходит полный pipeline.

Среда: свежая Release-сборка `build_olesnitskiy`, `g++-14`, `-O3 -DNDEBUG`, `std=gnu++23`; эти флаги видны в
[`build_olesnitskiy/compile_commands.json`](../../build_olesnitskiy/compile_commands.json). Число повторов задается framework-ом
`Perf::TaskRun`; отчет опирается на напечатанные строки `backend:task_run:time`.

## Сводка корректности

Источник функциональных тестов регистрирует все пять backend-ов
([`tests/functional/main.cpp`](tests/functional/main.cpp#L80)). Запуск
`PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests --gtest_filter='*olesnitskiy_v_hoare_sort_simple_merge*'`
выполнил 75 тестов: 60 passed для `seq/omp/stl/tbb`, 15 ALL skipped вне MPI. Отдельный запуск
`mpirun -np 2 ./build_olesnitskiy/bin/ppc_func_tests --gtest_filter='*olesnitskiy_v_hoare_sort_simple_merge_all_enabled*'`
прошел 15/15 ALL-тестов.

## Агрегированные результаты

| backend | time | speedup | efficiency | notes |
|---|---:|---:|---:|---|
| seq | 0.0058254364 s | 1.000 | 1.000 | `TaskRun`, `N=100000`, baseline |
| omp, 1 thread | 0.2042550788 s | 0.029 | 0.029 | много мелких блоков `kBlockSize=64` |
| omp, 2 threads | 0.2022266690 s | 0.029 | 0.014 | `PPC_NUM_THREADS=2` |
| omp, 4 threads | 0.1880041056 s | 0.031 | 0.008 | лучший OMP-замер |
| tbb, 1 worker | 0.0024417256 s | 2.386 | 2.386 | `blocked_range`, `parallel_for` |
| tbb, 2 workers | 0.0014976288 s | 3.889 | 1.945 | `PPC_NUM_THREADS=2` |
| tbb, 4 workers | 0.0011774682 s | 4.947 | 1.237 | лучший потоковый backend |
| stl, auto workers (12 на тестовой машине) | 0.0047718214 s | 1.221 | 0.102 | `hardware_concurrency`, `PPC_NUM_THREADS` не используется STL-кодом |
| all, 1 rank x 12 threads | 0.0126647332 s | 0.460 | 0.038 | total_workers=12 |
| all, 2 ranks x 12 threads | 0.0088037862 s | 0.662 | 0.028 | total_workers=24 |
| all, 4 ranks x 12 threads | 0.0043261284 s | 1.347 | 0.028 | total_workers=48 |

Дополнительные pipeline-замеры:

| backend | time | speedup | efficiency | notes |
|---|---:|---:|---:|---|
| seq | 0.0068995056 s | 1.000 | 1.000 | `pipeline`, baseline |
| omp, 1 thread | 0.1145523454 s | 0.060 | 0.060 | `pipeline` |
| omp, 4 threads | 0.1193738506 s | 0.058 | 0.014 | `pipeline` |
| stl, auto workers (12 на тестовой машине) | 0.0055675704 s | 1.239 | 0.103 | `pipeline`, `PPC_NUM_THREADS=1`, фактически auto-workers |
| tbb, 1 worker | 0.0067520042 s | 1.022 | 1.022 | `pipeline` |
| tbb, 4 workers | 0.0026312252 s | 2.622 | 0.656 | `pipeline` |
| all, 4 ranks x 12 threads | 0.0516886980 s | 0.133 | 0.003 | `pipeline`, total_workers=48 |

## Интерпретация различий

`seq` сортирует весь массив одним quicksort и служит baseline ([`seq/src/ops_seq.cpp`](seq/src/ops_seq.cpp#L75)).
`omp` создает параллельные области для множества блоков по 64 и выделяет временные векторы при слиянии
([`omp/src/ops_omp.cpp`](omp/src/ops_omp.cpp#L115), [`omp/src/ops_omp.cpp`](omp/src/ops_omp.cpp#L136)); на `N=100000`
измеренный speedup меньше 1. `tbb` использует `parallel_for(blocked_range)` для блоков и слияний
([`tbb/src/ops_tbb.cpp`](tbb/src/ops_tbb.cpp#L117)); при 4 workers получено `4.947x`. `stl` вручную создает и join-ит
потоки ([`stl/src/ops_stl.cpp`](stl/src/ops_stl.cpp#L38)); при 12 workers получено `1.221x`, но efficiency `0.102`.
`all` добавляет `MPI_Scatterv/Gatherv/Bcast` и финальное слияние на rank 0
([`all/src/ops_all.cpp`](all/src/ops_all.cpp#L246)); на 4 rank-ах время `0.0043261284 s`, но efficiency `0.028`.

## Репродуцируемость

Сборка:

```bash
cmake -S . -B build_olesnitskiy -DCMAKE_BUILD_TYPE=Release
cmake --build build_olesnitskiy --target ppc_func_tests ppc_perf_tests -j 4
```

Запуск корректности:

```bash
PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests --gtest_filter='*olesnitskiy_v_hoare_sort_simple_merge*'
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 mpirun -np 2 ./build_olesnitskiy/bin/ppc_func_tests --gtest_filter='*olesnitskiy_v_hoare_sort_simple_merge_all_enabled*'
```

Запуск замеров:

```bash
PPC_NUM_THREADS=1 ./build_olesnitskiy/bin/ppc_perf_tests --gtest_filter='*task_run*olesnitskiy_v_hoare_sort_simple_merge_seq_enabled'
PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_perf_tests --gtest_filter='*task_run*olesnitskiy_v_hoare_sort_simple_merge_tbb_enabled'
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 PPC_NUM_THREADS=1 mpirun -np 4 ./build_olesnitskiy/bin/ppc_perf_tests --gtest_filter='*task_run*olesnitskiy_v_hoare_sort_simple_merge_all_enabled'
PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_perf_tests --gtest_filter='*pipeline*olesnitskiy_v_hoare_sort_simple_merge_tbb_enabled'
```

## Заключение

Для измеренного размера `N=100000` лучшая численная версия — TBB с 4 workers: `0.0011774682 s`, speedup `4.947`.
ALL на 4 rank-ах быстрее baseline (`0.0043261284 s`, speedup `1.347`), но его efficiency низкая из-за 48 суммарных
workers и MPI-обменов. OMP на этом входе проигрывает baseline из-за накладных расходов.

## Источники

1. OpenMP API Specification 5.2, разделы data-sharing и schedule: https://www.openmp.org/spec-html/5.2/openmp.html
2. oneTBB `parallel_for`: https://www.intel.com/content/www/us/en/docs/onetbb/developer-guide-api-reference/2021-9/parallel-for.html
3. oneTBB `blocked_range`: https://oneapi-spec.uxlfoundation.org/specifications/oneapi/latest/elements/onetbb/source/algorithms/blocked_ranges/blocked_range_cls
4. MPI Standard, коллективные операции: https://www.mpi-forum.org/docs/
5. `std::thread::join`: https://en.cppreference.com/w/cpp/thread/thread/join.html
