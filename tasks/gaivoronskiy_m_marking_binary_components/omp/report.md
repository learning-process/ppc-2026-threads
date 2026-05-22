# Маркировка компонент на бинарном изображении — OMP

- Student: Гайворонский Максим Витальевич, group 3823Б1ПР1
- Technology: OMP
- Variant: 1

## 1. Контекст

OpenMP-версия распараллеливает strip-пайплайн маркировки компонент на бинарном изображении. Последовательная
реализация (`seq/report.md`) — эталон корректности и baseline для расчёта ускорения.

## 2. Постановка задачи

Входной вектор: `[rows, cols, pixels...]`.

- `0` — пиксель объекта;
- `1` — фон.

Выходной вектор того же размера: фон = `0`, каждая 4-связная компонента — уникальная положительная метка.

Ограничения проверяются в `ValidationImpl`: `rows > 0`, `cols > 0`, размер вектора равен `rows * cols + 2`.

## 3. Базовый алгоритм

Последовательная версия (`seq/src/ops_seq.cpp`) обходит изображение построчно и для каждого непромаркированного
пикселя объекта запускает BFS. Сложность по времени — O(rows × cols); по памяти — O(rows × cols) для выхода
и очереди BFS.

## 4. Схема распараллеливания

Число потоков: `num_threads = min(PPC_NUM_THREADS, rows)`. Число потоков задаётся через `PPC_NUM_THREADS`
(экспортируется также как `OMP_NUM_THREADS`).

| Этап                   | Директива OpenMP                                    | shared / private                                                                                | schedule               | Барьер                     |
| ---------------------- | --------------------------------------------------- | ----------------------------------------------------------------------------------------------- | ---------------------- | -------------------------- |
| Маркировка полос       | `#pragma omp parallel`                              | shared: input, output, labels_used, row_starts, cols; private: tid, next_label, локальные циклы | —                      | неявный в конце `parallel` |
| Смещение меток         | `#pragma omp parallel for`                          | shared: output, row_starts, base                                                                | `static`               | неявный                    |
| Union-find на границах | последовательно                                     | —                                                                                               | —                      | —                          |
| Первые позиции меток   | `#pragma omp parallel` + `#pragma omp parallel for` | shared: thread_first, first_pos; private: tid, local                                            | `static` (второй цикл) | неявный                    |
| Финальный remap        | `#pragma omp parallel for`                          | shared: input, output, parent, remap                                                            | `static`               | неявный                    |

**`default(none)`** используется во всех регионах — все переменные явно перечислены в `shared`.

**SEQ fallback:** при `num_threads <= 1` вызывается `RunSequentialMarking` без strip-пайплайна.

Ключевой фрагмент — параллельная маркировка полос:

```cpp
// File: omp/src/ops_omp.cpp
#pragma omp parallel num_threads(num_threads) default(none) \
    shared(input, output, labels_used, row_starts, cols, num_threads)
{
  const int tid = omp_get_thread_num();
  const int r_begin = row_starts[static_cast<std::size_t>(tid)];
  const int r_end = row_starts[static_cast<std::size_t>(tid) + 1];
  // BFS внутри полосы, запись сразу в output
  labels_used[static_cast<std::size_t>(tid)] = next_label;
}
```

Каждый поток пишет в непересекающиеся строки `output` — гонок записи нет. `labels_used[tid]` — независимые
ячейки, reduction не требуется.

**Выбор `schedule(static)`** для `parallel for` по полосам, корням и строкам: нагрузка предсказуема
(равномерные полосы / строки), статическое распределение минимизирует overhead планировщика.

## 5. Детали реализации

Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.

Pipeline `RunImpl`:

1. `ParallelLabelStrips` — BFS внутри полосы (`BfsLabelInStrip`), очередь на `std::vector`.
2. `AddBaseOffsets` — `#pragma omp parallel for schedule(static)`.
3. `MergeBoundariesUnionFind` + `FlattenParentForest` — последовательно.
4. `ComputeFirstPositionsParallel` — два OMP-региона (сбор локальных минимумов + merge по корням).
5. `BuildRemapFromFirstPositions` — последовательно.
6. `ApplyRemapInParallel` — `#pragma omp parallel for schedule(static)` по строкам, только пиксели объектов.

Оптимизации относительно базовой версии: маркировка сразу в `output` (без `local_planes`), `AddBaseOffsets`
вместо полного копирования полос, параллельный remap.

## 6. Проверка корректности

| Тест                   | Описание                     |
| ---------------------- | ---------------------------- |
| `single_L_component`   | L-образная компонента 3×3    |
| `four_separate_pixels` | четыре изолированных объекта |
| `all_background`       | только фон                   |
| `all_objects`          | только объекты               |
| `two_horizontal_bars`  | две горизонтальные полосы    |

Результат: **5/5** функциональных тестов пройдено. Сравнение с SEQ — побайтовое совпадение выхода.

## 7. Экспериментальная среда

| Параметр   | Значение                                          |
| ---------- | ------------------------------------------------- |
| CPU        | Apple M4, 10 ядер                                 |
| ОС         | macOS 15.5                                        |
| Компилятор | Apple Clang 17.0.0 + libomp                       |
| Сборка     | Release, `USE_FUNC_TESTS=ON`, `USE_PERF_TESTS=ON` |
| Perf-вход  | 2000×2000, чётные строки — объект                 |
| Повторы    | медиана по 3 сериям                               |
| Переменные | `PPC_NUM_THREADS` / `OMP_NUM_THREADS`             |

Команды:

```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
cmake --build build -j$(sysctl -n hw.ncpu)

./build/bin/ppc_func_tests --gtest_filter='*gaivoronskiy_m_marking_binary_components_omp*'

PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/task_run_gaivoronskiy_m_marking_binary_components_omp_enabled'
```

## 8. Результаты

Baseline SEQ — из `seq/report.md` (медиана 3 серий, perf-вход 2000×2000):

| Потоки   | SEQ, с   |
| -------- | -------- |
| 1        | 0.002901 |
| 2        | 0.002711 |
| 4        | 0.002683 |

Speedup = `T_seq / T_omp` при том же числе потоков (см. корневой `report.md`).

### task_run

| Потоки   | OMP, с   | Speedup vs SEQ   | Efficiency   |
| -------- | -------- | ---------------- | ------------ |
| 1        | 0.002615 | 1.11             | 111%         |
| 2        | 0.001753 | 1.55             | 77%          |
| 4        | 0.001388 | **1.93**         | **48%**      |

### pipeline (4 потока)

| Backend   | Время, с   | Speedup vs SEQ   |
| --------- | ---------- | ---------------- |
| SEQ       | 0.007558   | 1.00             |
| OMP       | 0.003381   | 2.24             |

На 500×500 OMP не обгонял SEQ из‑за overhead параллельных регионов. На 2000×2000 ускорение **~1.9×** при
4 потоках; `pipeline` — **~2.2×**.

## 9. Выводы

OpenMP даёт компактную реализацию strip-маркировки: `#pragma omp parallel` для полос, `parallel for schedule(static)`
для смещений и remap. `default(none)` заставляет явно описывать области видимости переменных. На больших
изображениях оптимизированная OMP-версия устойчиво быстрее SEQ; дальнейший рост потоков ограничен
последовательными этапами union-find и merge remap.
