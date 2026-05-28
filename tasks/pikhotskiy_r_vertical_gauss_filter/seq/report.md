# Линейная фильтрация (вертикальное разбиение), Гаусс 3x3 — SEQ

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: SEQ
- Variant: 25

## 1. Introduction

SEQ-версия реализована как baseline для сравнения с OMP, TBB, STL и ALL.

## 2. Problem Statement

- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`).
- Условие корректности входа: `data.size() == width * height`.
- Ограничения: `width > 0`, `height > 0`.
- Выход: изображение той же размерности после фильтра Гаусса 3x3.

## 3. Baseline Algorithm (Sequential)

Используется сепарабельная свертка:

1. Вертикальный проход `[1, 2, 1]` в промежуточный `int`-буфер.
2. Горизонтальный проход `[1, 2, 1]` в итоговый `uint8_t`-буфер.
3. Нормализация: `(sum + 15) / 16`.
4. Границы: `clamp` индексов.

## 4. Parallelization Scheme

Параллелизм отсутствует.
Обработка выполняется последовательно по вертикальным полосам.
Это и есть честный baseline.

## 5. Implementation Details

- Файлы: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterSEQ`.
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
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

Локально запускались функциональные и performance-тесты.

## 7. Results and Discussion

### 7.1 Correctness

Корректность подтверждается функциональными тестами.
Отдельно проверяются негативные случаи валидации входа.

### 7.2 Performance

Единые определения метрик:

- `workers` — число исполнительных единиц (для SEQ всегда `1`).
- `time` — wall-clock время, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 |   T_seq |    1.00 |           N/A |

## 8. Conclusions

SEQ-версия используется как эталон корректности.
Также это базовое время для расчета ускорения.

## 9. References

1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. Course report requirements: `docs/common_information/report.rst`
