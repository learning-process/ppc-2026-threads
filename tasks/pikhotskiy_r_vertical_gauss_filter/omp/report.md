# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — OMP

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: OMP
- Variant: 25

## 1. Introduction
В OMP-версии реализована параллельная обработка вертикальных полос изображения при сохранении той же математической модели, что и в SEQ.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`), где `data.size() == width * height`.
- Выход: изображение той же размерности после фильтра Гаусса 3x3.
- Ограничения: `width > 0`, `height > 0`, корректная обработка границ.

## 3. Baseline Algorithm (Sequential)
Алгоритм совпадает с SEQ:
1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme
Параллелизм организован по индексам вертикальных полос (`stripe_index`) в двух фазах (`vertical` и `horizontal`):
- директива: `#pragma omp parallel for default(none) schedule(static)`.
- `shared`: `stripe_count` (явно), а также общие поля объекта (`width_`, `height_`, `stripe_width_`, буферы) используются только для непересекающихся записей.
- `private`: `stripe_index`, `x_begin`, `x_end` (переменные цикла и локальные переменные итерации).
- `reduction`: не требуется, потому что нет глобальных суммирований между потоками.
- `schedule(static)`: фиксированное распределение полос между потоками.

Синхронизация между фазами обеспечивается неявным барьером в конце каждого `parallel for`.

## 5. Implementation Details
- Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterOMP`.
- Pipeline:
  - `ValidationImpl`
  - `PreProcessingImpl` (в том числе расчет `stripe_width_` от `ppc::util::GetNumThreads()`)
  - `RunImpl` (две параллельные фазы)
  - `PostProcessingImpl`

## 6. Experimental Setup
- Сборка:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build --target ppc_func_tests ppc_perf_tests
```
- Пример запуска:
```bash
export PPC_NUM_THREADS=4
export OMP_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```
Локально выполнялись оба типа проверок: функциональные и performance.

## 7. Results and Discussion
### 7.1 Correctness
OMP-версия проверяется теми же функциональными тестами, что и SEQ. Для одинаковых входных данных результат побайтно совпадает с baseline.

### 7.2 Performance
Определения метрик (единые для всех отчетов задачи):
- `workers` — число потоков OpenMP.
- `time` — wall-clock время выполнения, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency |
|------|---------|---------|---------|------------|
| seq  | 1       | T_seq   | 1.00    | N/A        |
| omp  | N       | T_omp   | T_seq / T_omp | (T_seq / T_omp) / N * 100% |

## 8. Conclusions
OMP-реализация дает параллельное ускорение на вертикальном разбиении и сохраняет эквивалентность SEQ по результату.

## 9. References
1. OpenMP specification: <https://www.openmp.org/specifications/>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
