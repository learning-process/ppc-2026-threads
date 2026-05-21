# Линейная фильтрация (вертикальное разбиение), Гаусс 3x3 — STL

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: STL
- Variant: 25

## 1. Introduction

STL-версия использует `std::thread` без внешнего runtime.

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

Изображение делится на вертикальные полосы.
Полосы распределяются между потоками.

Для реального параллелизма соблюден порядок:

1. Сначала запускаются все потоки текущей фазы.
2. Затем вызывается `join()` для всех запущенных потоков.
3. После завершения фазы запускается следующая фаза.

Это исключает последовательный шаблон "создал-подождал".

## 5. Implementation Details

- Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterSTL`.
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

STL-версия проходит те же функциональные проверки, что и baseline SEQ.

### 7.2 Performance

Единые определения метрик:

- `workers` — число потоков `std::thread`.
- `time` — wall-clock время, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 |   T_seq |    1.00 |           N/A |
| stl  |       N |   T_stl |   S_stl |         E_stl |

## 8. Conclusions

STL-версия дает переносимый контроль над потоками и корректное
параллельное выполнение при запуске всех потоков до `join()`.

## 9. References

1. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
