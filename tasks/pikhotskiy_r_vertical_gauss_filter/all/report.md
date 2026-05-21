# Линейная фильтрация (вертикальное разбиение), Гаусс 3x3 — ALL

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: ALL
- Variant: 25

## 1. Introduction

ALL-версия объединяет общий pipeline и единый backend запуска задачи.

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

Для `threads`-задачи используется конфигурация:
`ranks x threads = 1 x N`.

- `ranks = 1`: межпроцессного обмена нет.
- `threads = N`: обработка полос идет потоками.

Смысл MPI-синхронизации здесь:

- отдельные MPI-барьеры не применяются, так как процесс один;
- синхронизация обеспечивается границей фаз
  (`vertical` полностью завершается до `horizontal`).

## 5. Implementation Details

- Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterALL`.
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
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

Локально запускались функциональные и performance-тесты.

## 7. Results and Discussion

### 7.1 Correctness

ALL-версия проходит те же функциональные проверки, что и baseline SEQ.

### 7.2 Performance

Единые определения метрик:

- `workers` — число потоков в процессе (`workers = threads`).
- `time` — wall-clock время, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 |   T_seq |    1.00 |           N/A |
| all  |       N |   T_all |   S_all |         E_all |

## 8. Conclusions

ALL-версия поддерживает единый запуск в инфраструктуре курса
и сохраняет корректность baseline.

## 9. References

1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. Course report requirements: `docs/common_information/report.rst`
