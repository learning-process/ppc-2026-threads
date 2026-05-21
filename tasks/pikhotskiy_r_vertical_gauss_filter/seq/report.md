# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — SEQ

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: SEQ
- Variant: 25

## 1. Introduction
Цель реализации SEQ-версии: получить корректный baseline для фильтра Гаусса 3x3 с вертикальным разбиением, относительно которого сравниваются OMP/TBB/STL/ALL.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`), где `data.size() == width * height`.
- Выход: изображение той же размерности после фильтрации ядром Гаусса 3x3.
- Ограничения: `width > 0`, `height > 0`, корректная обработка границ.

## 3. Baseline Algorithm (Sequential)
Используется сепарабельная свертка:
1. Вертикальный проход с коэффициентами `[1, 2, 1]` в промежуточный `int`-буфер.
2. Горизонтальный проход с коэффициентами `[1, 2, 1]` в итоговый `uint8_t`-буфер.
3. Нормализация: `(sum + 15) / 16`.
4. Границы обрабатываются через `clamp` индексов.

## 4. Parallelization Scheme
Параллелизм не используется. Обработка идет последовательно по вертикальным полосам фиксированной ширины (`kStripeDivider = 8`), что задает структуру baseline и совпадает с логикой декомпозиции в параллельных версиях.

## 5. Implementation Details
- Файлы: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterSEQ`.
- Pipeline:
  - `ValidationImpl`: проверка размеров и `data.size()`.
  - `PreProcessingImpl`: подготовка `source_`, `vertical_buffer_`, `result_buffer_`.
  - `RunImpl`: вертикальный и горизонтальный проходы.
  - `PostProcessingImpl`: запись результата в `OutType`.

## 6. Experimental Setup
- Сборка:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build --target ppc_func_tests ppc_perf_tests
```
- Запуск функциональных тестов:
```bash
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```
- Запуск performance-тестов:
```bash
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```
Локально выполнялись оба типа проверок: функциональные и performance.

## 7. Results and Discussion
### 7.1 Correctness
Корректность проверяется функциональными тестами на валидных и невалидных входах, а также совпадением ожидаемых значений после фильтрации.

### 7.2 Performance
Определения метрик (единые для всех отчетов задачи):
- `workers` — число исполнительных единиц (для SEQ всегда `1`).
- `time` — wall-clock время выполнения, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency |
|------|---------|---------|---------|------------|
| seq  | 1       | T_seq   | 1.00    | N/A        |

## 8. Conclusions
SEQ-реализация является честным baseline без параллелизма и используется как эталон корректности и базовое время для расчета ускорения в параллельных версиях.

## 9. References
1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. Course report requirements: `docs/common_information/report.rst`
