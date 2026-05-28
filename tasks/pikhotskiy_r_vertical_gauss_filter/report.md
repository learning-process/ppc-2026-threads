# Линейная фильтрация (вертикальное разбиение), Гаусс 3x3

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Task: `pikhotskiy_r_vertical_gauss_filter`
- Variant: 25
- Implementations: SEQ, OMP, TBB, STL, ALL

## 1. Introduction

Реализован фильтр Гаусса 3x3 с вертикальным разбиением.
Поддержаны пять технологий:
`seq`, `omp`, `tbb`, `stl`, `all`.

## 2. Problem Statement

- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`).
- Условие корректности входа: `data.size() == width * height`.
- Ограничения: `width > 0`, `height > 0`.
- Выход: изображение после фильтра Гаусса 3x3.

## 3. Baseline Algorithm (Sequential)

Общий алгоритм для всех реализаций:

1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme

- `SEQ`: без параллелизма, baseline.
- `OMP`: `parallel for` по полосам, `schedule(static)`.
- `TBB`: `parallel_for` + `blocked_range`.
- `STL`: `std::thread` по полосам.
- `ALL`: единый backend; для threads-конфигурации `1 x N`.

## 5. Implementation Details

- Общие типы: `common/include/common.hpp`.
- Папки реализаций: `seq`, `omp`, `tbb`, `stl`, `all`.
- Единый pipeline: `ValidationImpl`, `PreProcessingImpl`, `RunImpl`,
  `PostProcessingImpl`.

## 6. Experimental Setup

Сборка:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build --target ppc_func_tests ppc_perf_tests
```

Запуск:

```bash
export PPC_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

Локально запускались функциональные и performance-тесты.

## 7. Results and Discussion

### 7.1 Correctness

Проверка корректности выполняется функциональными тестами.
Параллельные реализации сравниваются с baseline SEQ.

### 7.2 Performance

Единые определения метрик:

- `workers` — число исполнительных единиц.
- `time` — wall-clock время, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 |   T_seq |    1.00 |           N/A |
| omp  |       N |   T_omp |   S_omp |         E_omp |
| tbb  |       N |   T_tbb |   S_tbb |         E_tbb |
| stl  |       N |   T_stl |   S_stl |         E_stl |
| all  |       N |   T_all |   S_all |         E_all |

## 8. Conclusions

Задача реализована и документирована в едином стиле.
Это сделано для всех технологий.
SEQ используется как baseline, остальные версии сравниваются с ним.

## 9. References

1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. OpenMP specification: <https://www.openmp.org/specifications/>
3. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
4. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
5. Course report requirements: `docs/common_information/report.rst`
