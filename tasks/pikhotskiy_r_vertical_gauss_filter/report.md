# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3

- Student: Pikhotskiy Roman Vladimirovich, group 3823B1FI1
- Task: `pikhotskiy_r_vertical_gauss_filter`
- Variant: 25
- Implementations: SEQ, OMP, TBB, STL, ALL

## 1. Introduction
В задаче реализован фильтр Гаусса 3x3 для изображений с вертикальным разбиением и единым `BaseTask` pipeline.  
Подготовлены отчеты по всем требуемым технологиям: `seq`, `omp`, `tbb`, `stl`, `all`.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`std::vector<std::uint8_t>`), где `data.size() == width * height`.
- Выход: отфильтрованное изображение той же размерности.
- Обязательные проверки:
  - `width > 0`, `height > 0`;
  - соответствие размера данных;
  - корректная обработка границ.

## 3. Baseline Algorithm (Sequential)
Во всех версиях используется один и тот же базовый алгоритм:
1. Вертикальный проход ядром `[1, 2, 1]`.
2. Горизонтальный проход ядром `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme
- `SEQ`: baseline без параллелизма.
- `OMP`: `parallel for` по полосам, `schedule(static)`.
- `TBB`: `parallel_for` + `blocked_range`, управление задачами runtime oneTBB.
- `STL`: `std::thread` с явным запуском worker-потоков по полосам.
- `ALL`: универсальный backend для общей инфраструктуры; для `threads`-задачи конфигурация рассматривается как `1 x N` (один процесс, N потоков).

## 5. Implementation Details
- Общие типы: `common/include/common.hpp`.
- Реализации:
  - `seq/include|src`
  - `omp/include|src`
  - `tbb/include|src`
  - `stl/include|src`
  - `all/include|src`
- Единый pipeline в каждом backend:
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
- Запуск:
```bash
export PPC_NUM_THREADS=4
./build/bin/ppc_func_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```
Локально выполнялись оба типа проверок: функциональные и performance.

## 7. Results and Discussion
### 7.1 Correctness
Проверка корректности выполняется функциональными тестами, включая проверку валидации входа и сравнение результатов параллельных версий с baseline SEQ.

### 7.2 Performance
Единые определения метрик для всех отчетов задачи:
- `workers` — число исполнительных единиц (потоки; для SEQ: `1`).
- `time` — wall-clock время выполнения, секунды.
- `speedup = T_seq / T_mode`.
- `efficiency = speedup / workers * 100%`.

| mode | workers | time, s | speedup | efficiency |
|------|---------|---------|---------|------------|
| seq  | 1       | T_seq   | 1.00    | N/A        |
| omp  | N       | T_omp   | T_seq / T_omp | (T_seq / T_omp) / N * 100% |
| tbb  | N       | T_tbb   | T_seq / T_tbb | (T_seq / T_tbb) / N * 100% |
| stl  | N       | T_stl   | T_seq / T_stl | (T_seq / T_stl) / N * 100% |
| all  | N       | T_all   | T_seq / T_all | (T_seq / T_all) / N * 100% |

## 8. Conclusions
Задача оформлена в едином стиле для всех технологий, с согласованными определениями метрик и воспроизводимыми командами запуска.  
SEQ выступает baseline, остальные backend-ы сравниваются относительно него по корректности и производительности.

## 9. References
1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. OpenMP specification: <https://www.openmp.org/specifications/>
3. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
4. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
5. Course report requirements: `docs/common_information/report.rst`
