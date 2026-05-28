# Линейная фильтрация (вертикальное разбиение), Гаусс 3x3 — TBB

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: TBB
- Variant: 25

## 1. Introduction

TBB-версия использует task-based модель oneTBB для обработки вертикальных
полос изображения.

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

Используется `oneapi::tbb::parallel_for` по диапазону полос:
`blocked_range<int>(0, stripe_count, grainsize)`.

- `blocked_range`: задает блок индексов полос.
- `grainsize`: минимальный размер блока (не менее 1 полосы).
- `partitioner`: стандартный `auto_partitioner` у `parallel_for`.
- Контроль конкуренции: каждая задача пишет только в свою полосу.

Первая фаза (`vertical`) полностью завершается до запуска второй
(`horizontal`), поэтому между фазами нет гонок.

## 5. Implementation Details

- Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterTBB`.
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

TBB-версия проходит те же функциональные проверки, что и baseline SEQ.

### 7.2 Performance

Единые определения метрик:

- `workers` — число рабочих потоков TBB runtime.
- `time` — wall-clock время, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 |   T_seq |    1.00 |           N/A |
| tbb  |       N |   T_tbb |   S_tbb |         E_tbb |

## 8. Conclusions

TBB-версия сохраняет корректность baseline.
Распределение полос остается гибким за счет task-based исполнения.

## 9. References

1. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
