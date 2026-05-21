# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — TBB

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: TBB
- Variant: 25

## 1. Introduction
TBB-версия использует модель задач oneTBB для параллельной обработки вертикальных полос изображения.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`), где `data.size() == width * height`.
- Выход: изображение после фильтра Гаусса 3x3, размер сохраняется.
- Ограничения: `width > 0`, `height > 0`, корректная обработка границ.

## 3. Baseline Algorithm (Sequential)
Вычислительное ядро совпадает с SEQ:
1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme
Распараллеливание выполняется через `oneapi::tbb::parallel_for` по диапазону вертикальных полос:
- `blocked_range<int>(0, stripe_count, grainsize)`.
- `grainsize` выбирается как минимум `1` полоса; дальше oneTBB может объединять/делить задачи динамически.
- `partitioner`: используется стандартная стратегия oneTBB (`auto_partitioner` по умолчанию для `parallel_for`), чтобы балансировать загрузку.
- Контроль конкуренции: каждая задача пишет только в свою полосу буфера, пересечений записи нет.

Между вертикальной и горизонтальной фазами есть явная фазовая граница (сначала завершается весь первый `parallel_for`, затем стартует второй).

## 5. Implementation Details
- Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterTBB`.
- Pipeline:
  - `ValidationImpl`
  - `PreProcessingImpl`
  - `RunImpl`
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
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```
Локально выполнялись оба типа проверок: функциональные и performance.

## 7. Results and Discussion
### 7.1 Correctness
Корректность проверяется теми же функциональными тестами, что и для SEQ/OMP, с проверкой совпадения результата с baseline.

### 7.2 Performance
Определения метрик (единые для всех отчетов задачи):
- `workers` — число рабочих потоков, доступных TBB runtime.
- `time` — wall-clock время выполнения, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency |
|------|---------|---------|---------|------------|
| seq  | 1       | T_seq   | 1.00    | N/A        |
| tbb  | N       | T_tbb   | T_seq / T_tbb | (T_seq / T_tbb) / N * 100% |

## 8. Conclusions
TBB-версия сохраняет корректность baseline и использует более гибкое task-based распределение полос, что обычно улучшает утилизацию потоков при росте размера изображения.

## 9. References
1. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
