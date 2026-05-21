# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — STL

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Technology: STL
- Variant: 25

## 1. Introduction
STL-версия реализует параллелизм средствами `std::thread` без внешнего runtime.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`), где `data.size() == width * height`.
- Выход: изображение после фильтра Гаусса 3x3 с сохранением размеров.
- Ограничения: `width > 0`, `height > 0`, корректная обработка границ.

## 3. Baseline Algorithm (Sequential)
Алгоритм совпадает с SEQ:
1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme
Изображение делится на вертикальные полосы, полосы распределяются между потоками.

Критично для реального параллелизма:
1. Сначала запускаются все потоки для текущей фазы.
2. Только после запуска всех потоков вызывается `join()` для каждого потока.
3. После завершения фазы запускается следующая фаза (vertical -> horizontal).

Такой порядок исключает последовательный запуск "создал-подождал" и дает фактическое конкурентное выполнение.

## 5. Implementation Details
- Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.
- Класс: `PikhotskiyRVerticalGaussFilterSTL`.
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
STL-версия проходит те же функциональные проверки, что baseline, и должна совпадать с ним по выходным данным.

### 7.2 Performance
Определения метрик (единые для всех отчетов задачи):
- `workers` — число запущенных `std::thread`.
- `time` — wall-clock время выполнения, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency |
|------|---------|---------|---------|------------|
| seq  | 1       | T_seq   | 1.00    | N/A        |
| stl  | N       | T_stl   | T_seq / T_stl | (T_seq / T_stl) / N * 100% |

## 8. Conclusions
STL-версия дает переносимый контроль над потоками и корректный параллелизм при условии запуска всех worker-потоков до этапа `join()`.

## 9. References
1. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
