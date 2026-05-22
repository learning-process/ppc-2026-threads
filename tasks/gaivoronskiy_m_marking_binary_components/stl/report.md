# Маркировка компонент на бинарном изображении — STL

- Student: Гайворонский Максим Витальевич, group 3823Б1ПР1
- Technology: STL (`std::thread`)
- Variant: 1

## 1. Контекст

STL-версия реализует strip-пайплайн маркировки компонент с явным управлением потоками через `std::thread`.
Последовательная реализация (`seq/report.md`) — эталон корректности и baseline для ускорения.

## 2. Постановка задачи

Входной вектор: `[rows, cols, pixels...]`.

- `0` — пиксель объекта;
- `1` — фон.

Выходной вектор того же размера: фон = `0`, каждая 4-связная компонента — уникальная положительная метка.

Ограничения проверяются в `ValidationImpl`: `rows > 0`, `cols > 0`, размер вектора равен `rows * cols + 2`.

## 3. Базовый алгоритм

Последовательная версия (`seq/src/ops_seq.cpp`) — построчный обход и BFS для каждой новой компоненты.
Сложность по времени — O(rows × cols); по памяти — O(rows × cols).

## 4. Схема распараллеливания

Число потоков: `num_threads = min(PPC_NUM_THREADS, rows)`.

| Этап                 | Разбиение                                      | Локальные данные      | Синхронизация                               |
| -------------------- | ---------------------------------------------- | --------------------- | ------------------------------------------- |
| Маркировка полос     | `tid` → `[row_starts[tid], row_starts[tid+1])` | `labels_used[tid]`    | `join` после создания всех потоков          |
| Смещение меток       | по `thread_idx`                                | —                     | `join`                                      |
| Union-find           | последовательно                                | `parent[]`            | —                                           |
| Первые позиции меток | по `tid` → диапазон строк                      | `thread_first[tid][]` | `join`, затем последовательный merge корней |
| Финальный remap      | по `tid` → диапазон строк                      | —                     | `join`                                      |

**Паттерн create → work → join:** потоки `1 … num_threads-1` создаются через `emplace_back`, поток `0`
выполняет работу на главном потоке, затем все worker-ы `join()`. Синхронизация через `mutex`/`atomic`
не требуется: каждый поток пишет в непересекающиеся строки `output` и в свою ячейку `labels_used[tid]`.

**SEQ fallback:** при `num_threads <= 1` — `RunSequentialMarking` без strip-пайплайна.

Ключевой фрагмент:

```cpp
// File: stl/src/ops_stl.cpp
std::vector<std::thread> workers;
workers.reserve(static_cast<std::size_t>(std::max(0, num_threads - 1)));

auto worker = [&](int tid) {
  const int r_begin = row_starts[static_cast<std::size_t>(tid)];
  const int r_end = row_starts[static_cast<std::size_t>(tid) + 1];
  labels_used[static_cast<std::size_t>(tid)] = LabelSingleStrip(input, output, cols, r_begin, r_end);
};

for (int tid = 1; tid < num_threads; ++tid) {
  workers.emplace_back(worker, tid);
}
worker(0);
for (auto &thread : workers) {
  thread.join();
}
```

`join` вызывается **после** запуска всех потоков — реальный параллелизм сохраняется.

## 5. Детали реализации

Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.

Pipeline `RunImpl`:

1. `ParallelLabelStrips` — BFS в полосе (`BfsLabelInStrip`), очередь на `std::vector`.
2. `AddBaseOffsets` — параллельное добавление глобальных смещений.
3. `MergeBoundariesUnionFind` + `FlattenParentForest` — последовательно.
4. `ComputeFirstPositionsParallel` + `BuildRemapFromFirstPositions`.
5. `ApplyRemapParallel` — remap только пикселей объектов.

Оптимизации: маркировка сразу в `output` (без `local_planes`), `AddBaseOffsets` вместо копирования полос,
параллельный сбор first positions, SEQ fallback.

## 6. Проверка корректности

| Тест                   | Описание                     |
| ---------------------- | ---------------------------- |
| `single_L_component`   | L-образная компонента 3×3    |
| `four_separate_pixels` | четыре изолированных объекта |
| `all_background`       | только фон                   |
| `all_objects`          | только объекты               |
| `two_horizontal_bars`  | две горизонтальные полосы    |

Результат: **5/5** функциональных тестов пройдено. Выход совпадает с SEQ.

Гонок записи нет: потоки работают с непересекающимися диапазонами строк и независимыми `labels_used[tid]`.

## 7. Экспериментальная среда

| Параметр   | Значение                          |
| ---------- | --------------------------------- |
| CPU        | Apple M4, 10 ядер                 |
| ОС         | macOS 15.5                        |
| Компилятор | Apple Clang 17.0.0, Release       |
| Perf-вход  | 2000×2000, чётные строки — объект |
| Повторы    | медиана по 3 сериям               |
| Переменные | `PPC_NUM_THREADS`                 |

Команды:

```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
cmake --build build -j$(sysctl -n hw.ncpu)

./build/bin/ppc_func_tests --gtest_filter='*gaivoronskiy_m_marking_binary_components_stl*'

PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/task_run_gaivoronskiy_m_marking_binary_components_stl_enabled'
```

## 8. Результаты

Baseline SEQ — из `seq/report.md` (медиана 3 серий, perf-вход 2000×2000):

| Потоки   | SEQ, с   |
| -------- | -------- |
| 1        | 0.002901 |
| 2        | 0.002711 |
| 4        | 0.002683 |

Speedup = `T_seq / T_stl` при том же числе потоков (см. корневой `report.md`).

### task_run

| Потоки   | STL, с   | Speedup vs SEQ   | Efficiency   |
| -------- | -------- | ---------------- | ------------ |
| 1        | 0.002580 | 1.12             | 112%         |
| 2        | 0.001905 | 1.42             | 71%          |
| 4        | 0.001024 | **2.62**         | **66%**      |

### pipeline (4 потока)

| Backend   | Время, с   | Speedup vs SEQ   |
| --------- | ---------- | ---------------- |
| SEQ       | 0.007558   | 1.00             |
| STL       | 0.003651   | 2.07             |

На 500×500 STL проигрывал SEQ из‑за overhead `create/join`. На 2000×2000 speedup **~2.6×** при 4 потоках —
сопоставимо с TBB и выше OMP по `task_run`.

## 9. Выводы

Ручная реализация на `std::thread` требует явного разбиения данных и корректного `join`, но даёт полный
контроль над границами работы. При оптимизации памяти (маркировка в `output`, параллельный remap) STL
на больших изображениях показывает ускорение, сравнимое с TBB. Основной overhead — создание/уничтожение
потоков на каждом вызове `RunImpl`.
