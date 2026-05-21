# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — ALL

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: ALL
- Variant: 25

## 1. Introduction
ALL-версия объединяет общий pipeline задачи и предоставляет универсальный backend для запуска в общей инфраструктуре курса.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`), где `data.size() == width * height`.
- Выход: изображение после фильтра Гаусса 3x3 с сохранением размеров.
- Ограничения: `width > 0`, `height > 0`, корректная обработка границ.

## 3. Baseline Algorithm (Sequential)
Вычислительное ядро одинаково с SEQ:
1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme
Для задачи типа `threads` конфигурация выполнения интерпретируется как `ranks x threads = 1 x N`:
- `ranks = 1`, межпроцессного обмена нет.
- `threads = N`, обработка выполняется по вертикальным полосам.

Смысл MPI-синхронизации в этой конфигурации:
- отдельные MPI-барьеры не используются, так как расчет идет в одном процессе;
- роль синхронизации выполняют фазовые границы между проходами (после завершения vertical-pass запускается horizontal-pass).

## 5. Implementation Details
- Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterALL`.
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
ALL-версия проверяется теми же функциональными тестами и должна совпадать с baseline по выходным данным.

### 7.2 Performance
Определения метрик (единые для всех отчетов задачи):
- `workers` — число потоков в процессе (`workers = threads`, так как `ranks = 1`).
- `time` — wall-clock время выполнения, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency |
|------|---------|---------|---------|------------|
| seq  | 1       | T_seq   | 1.00    | N/A        |
| all  | N       | T_all   | T_seq / T_all | (T_seq / T_all) / N * 100% |

## 8. Conclusions
ALL-версия поддерживает единый запуск в framework-е курса, при этом сохраняет те же вычисления и проверяемую корректность, что и остальные реализации задачи.

## 9. References
1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. Course report requirements: `docs/common_information/report.rst`
