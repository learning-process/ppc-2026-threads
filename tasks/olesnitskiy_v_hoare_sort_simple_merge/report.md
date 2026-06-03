# Сортировка Хоара с простым слиянием — сводный отчет

- **Студент:** Олесницкий Владимир Тарасович, 3823Б1ПР2
- **Технологии:** SEQ, OMP, TBB, STL, ALL (MPI + STL)
- **Вариант:** 13

## 1. Контекст

Работа сравнивает реализации сортировки `std::vector<int>` по неубыванию для
backend-ов `seq`, `omp`, `tbb`, `stl` и `all`. Последовательная версия задает
baseline, потоковые версии проверяют разные модели внутрипроцессного
распараллеливания, а ALL добавляет MPI-уровень.

## 2. Постановка задачи

- **Входные данные:** непустой объект `std::vector<int>`.
- **Выходные данные:** объект `std::vector<int>` того же размера, элементы
  которого переставлены в порядке неубывания.
- **Эталон корректности:** совпадение с результатом `std::ranges::sort`.
- **Baseline:** `seq`, `PPC_NUM_THREADS=1`, `T_seq = 0.0058254364 s`.

Единые типы входа и выхода заданы в `common/include/common.hpp`.

## 3. Базовый алгоритм

Последовательное ядро использует сортировку Хоара. Массив разбивается по
опорному элементу из середины диапазона, затем поддиапазоны сортируются
итеративным quicksort.

Параллельные реализации используют общий принцип простого слияния: сначала
сортируются независимые блоки по 64 элемента, затем уровни слияния объединяют
соседние отсортированные диапазоны.

## 4. Схемы распараллеливания

- **OMP:** два `parallel for` для сортировки блоков и слияния пар диапазонов.
- **TBB:** `oneapi::tbb::parallel_for` с `blocked_range<size_t>`.
- **STL:** ручной запуск `std::thread`; число потоков берется из
  `std::thread::hardware_concurrency()`.
- **ALL:** MPI распределяет фрагменты между rank-ами, а локальная часть каждого
  rank-а использует STL-сортировку.

Для OMP/TBB применялись `PPC_NUM_THREADS=1,2,4`; TBB дополнительно
ограничивается `tbb::global_control`. STL и локальная часть ALL берут число
потоков из `std::thread::hardware_concurrency()`. На тестовой машине `nproc =
12`, поэтому такие строки обозначены как auto workers.

## 5. Детали методики

Performance-тест генерирует `N=100000` случайных `int` из диапазона
`[-1000000, 1000000]`. Использован `TaskRun`; режим выбирается в раннере.
Speedup считался как `T_seq / T_backend`, efficiency — как
`speedup / workers`.

Ограничение методики: `TaskRun` повторяет только `Run()` после одного
`PreProcessing()`, а вход создается заново для каждого gtest-параметра через
`std::random_device`. Поэтому дополнительно сняты pipeline-замеры, где каждый
повтор проходит полный pipeline.

## 6. Проверка корректности

Функциональные тесты регистрируют все пять backend-ов. Запуск
`PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests` с фильтром
`*olesnitskiy_v_hoare_sort_simple_merge*` выполнил 75 тестов: 60 passed для
`seq/omp/stl/tbb`, 15 ALL skipped вне MPI.

Отдельный запуск `mpirun -np 2 ./build_olesnitskiy/bin/ppc_func_tests` с
фильтром `*olesnitskiy_v_hoare_sort_simple_merge_all_enabled*` прошел 15/15
ALL-тестов.

## 7. Экспериментальная среда

- **Сборка:** `build_olesnitskiy`
- **Compiler:** `g++-14`
- **Flags:** `-O3 -DNDEBUG`, `std=gnu++23`
- **Размер входных данных:** `N=100000`
- **Диапазон значений:** `[-1000000, 1000000]`
- **Число повторов:** 5 по умолчанию
- **Auto workers:** 12 на тестовой машине

## 8. Результаты TaskRun

- backend: seq; time: 0.0058254364 s; speedup: 1.000; efficiency: 1.000; notes:
  `TaskRun`, `N=100000`, baseline.
- backend: omp, 1 thread; time: 0.2042550788 s; speedup: 0.029; efficiency:
  0.029; notes: много мелких блоков `kBlockSize=64`.
- backend: omp, 2 threads; time: 0.2022266690 s; speedup: 0.029; efficiency:
  0.014; notes: `PPC_NUM_THREADS=2`.
- backend: omp, 4 threads; time: 0.1880041056 s; speedup: 0.031; efficiency:
  0.008; notes: лучший OMP-замер.
- backend: tbb, 1 worker; time: 0.0024417256 s; speedup: 2.386; efficiency:
  2.386; notes: `blocked_range`, `parallel_for`.
- backend: tbb, 2 workers; time: 0.0014976288 s; speedup: 3.889; efficiency:
  1.945; notes: `PPC_NUM_THREADS=2`.
- backend: tbb, 4 workers; time: 0.0011774682 s; speedup: 4.947; efficiency:
  1.237; notes: лучший потоковый backend.
- backend: stl, auto workers (12 на тестовой машине); time: 0.0047718214 s;
  speedup: 1.221; efficiency: 0.102; notes: `hardware_concurrency`; без env.
- backend: all, 1 rank x 12 threads; time: 0.0126647332 s; speedup: 0.460;
  efficiency: 0.038; notes: total_workers=12.
- backend: all, 2 ranks x 12 threads; time: 0.0088037862 s; speedup: 0.662;
  efficiency: 0.028; notes: total_workers=24.
- backend: all, 4 ranks x 12 threads; time: 0.0043261284 s; speedup: 1.347;
  efficiency: 0.028; notes: total_workers=48.

## 9. Результаты PipelineRun

- backend: seq; time: 0.0068995056 s; speedup: 1.000; efficiency: 1.000; notes:
  `pipeline`, baseline.
- backend: omp, 1 thread; time: 0.1145523454 s; speedup: 0.060; efficiency:
  0.060; notes: `pipeline`.
- backend: omp, 4 threads; time: 0.1193738506 s; speedup: 0.058; efficiency:
  0.014; notes: `pipeline`.
- backend: stl, auto workers (12 на тестовой машине); time: 0.0055675704 s;
  speedup: 1.239; efficiency: 0.103; notes: `pipeline`; auto.
- backend: tbb, 1 worker; time: 0.0067520042 s; speedup: 1.022; efficiency:
  1.022; notes: `pipeline`.
- backend: tbb, 4 workers; time: 0.0026312252 s; speedup: 2.622; efficiency:
  0.656; notes: `pipeline`.
- backend: all, 4 ranks x 12 threads; time: 0.0516886980 s; speedup: 0.133;
  efficiency: 0.003; notes: `pipeline`, total_workers=48.

## 10. Интерпретация результатов

`seq` сортирует весь массив одним quicksort и служит baseline. `omp` создает
параллельные области для множества блоков по 64 и выделяет временные векторы при
слиянии; на `N=100000` измеренный speedup меньше 1.

`tbb` использует `parallel_for(blocked_range)` для блоков и слияний; при 4
workers получено `4.947x`. `stl` вручную создает и join-ит потоки; при 12
workers получено `1.221x`, но efficiency равна `0.102`.

`all` добавляет `MPI_Scatterv/Gatherv/Bcast` и финальное слияние на rank 0. На
4 rank-ах время `0.0043261284 s`, но efficiency равна `0.028`.

## 11. Репродуцируемость

Сборка:

```bash
cmake -S . -B build_olesnitskiy -DCMAKE_BUILD_TYPE=Release
cmake --build build_olesnitskiy --target ppc_func_tests ppc_perf_tests -j 4
```

Запуск корректности:

```bash
PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_func_tests \
  --gtest_filter='*olesnitskiy_v_hoare_sort_simple_merge*'
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  mpirun -np 2 ./build_olesnitskiy/bin/ppc_func_tests \
  --gtest_filter='*olesnitskiy_v_hoare_sort_simple_merge_all_enabled*'
```

Запуск замеров:

```bash
PPC_NUM_THREADS=1 ./build_olesnitskiy/bin/ppc_perf_tests \
  --gtest_filter='*task_run*olesnitskiy_v_hoare_sort_simple_merge_seq_enabled'
PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_perf_tests \
  --gtest_filter='*task_run*olesnitskiy_v_hoare_sort_simple_merge_tbb_enabled'
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 PPC_NUM_THREADS=1 \
  mpirun -np 4 ./build_olesnitskiy/bin/ppc_perf_tests \
  --gtest_filter='*task_run*olesnitskiy_v_hoare_sort_simple_merge_all_enabled'
PPC_NUM_THREADS=4 ./build_olesnitskiy/bin/ppc_perf_tests \
  --gtest_filter='*pipeline*olesnitskiy_v_hoare_sort_simple_merge_tbb_enabled'
```

## 12. Выводы

Для измеренного размера `N=100000` лучшая численная версия — TBB с 4 workers:
`0.0011774682 s`, speedup `4.947`. ALL на 4 rank-ах быстрее baseline
(`0.0043261284 s`, speedup `1.347`), но его efficiency низкая из-за 48
суммарных workers и MPI-обменов. OMP на этом входе проигрывает baseline из-за
накладных расходов.
