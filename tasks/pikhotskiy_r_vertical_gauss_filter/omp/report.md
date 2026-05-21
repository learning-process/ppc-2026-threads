# Линейная фильтрация (вертикальное разбиение), Гаусс 3x3 — OMP

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: OMP
- Variant: 25

## 1. Introduction

OMP-версия повторяет математику SEQ.
Распараллеливание идет по вертикальным полосам.

## 2. Problem Statement

- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`).
- Условие корректности входа: `data.size() == width * height`.
- Ограничения: `width > 0`, `height > 0`.
- Выход: изображение после фильтра Гаусса 3x3.

## 3. Baseline Algorithm (Sequential)

Алгоритм совпадает с SEQ:

1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme

Используется `#pragma omp parallel for default(none) schedule(static)` для
итераций по полосам.

- `shared`: `stripe_count` и общие буферы объекта (без пересечения записи).
- `private`: `stripe_index`, `x_begin`, `x_end`.
- `reduction`: не используется, глобального суммирования нет.
- `schedule(static)`: фиксированное распределение полос.

После каждого `parallel for` действует неявный барьер, поэтому фазы
`vertical` и `horizontal` синхронизированы.

## 5. Implementation Details

- Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterOMP`.
- Pipeline: `ValidationImpl`, `PreProcessingImpl`, `RunImpl`,
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
export OMP_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

Локально запускались функциональные и performance-тесты.

## 7. Results and Discussion

### 7.1 Correctness

OMP-версия проходит те же функциональные проверки, что и baseline SEQ.

### 7.2 Performance

Единые определения метрик:

- `workers` — число потоков OpenMP.
- `time` — wall-clock время, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 |   T_seq |    1.00 |           N/A |
| omp  |       N |   T_omp |   S_omp |         E_omp |

## 8. Conclusions

OMP-версия сохраняет корректность baseline и обеспечивает ускорение
за счет распараллеливания по вертикальным полосам.

## 9. References

1. OpenMP specification: <https://www.openmp.org/specifications/>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
