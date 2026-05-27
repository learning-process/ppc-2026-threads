# Маркировка компонент на бинарном изображении — TBB

- Student: Гайворонский Максим Витальевич, group 3823Б1ПР1
- Technology: TBB
- Variant: 1

## 1. Контекст

Задача — маркировка 4-связных компонент объектов на бинарном изображении. TBB-версия распараллеливает
strip-пайплайн через `tbb::parallel_for` и ограничивает конкуренцию через `tbb::task_arena`. Последовательная
реализация (`seq/report.md`) используется как эталон корректности и baseline для ускорения.

## 2. Постановка задачи

Входной вектор: `[rows, cols, pixels...]`.

- `0` — пиксель объекта;
- `1` — фон.

Выходной вектор того же размера: фон = `0`, каждая 4-связная компонента объекта — уникальная положительная метка.

Ограничения проверяются в `ValidationImpl`: `rows > 0`, `cols > 0`, размер вектора равен `rows * cols + 2`.

## 3. Базовый алгоритм

Последовательная версия (`seq/src/ops_seq.cpp`) обходит изображение построчно и для каждого непромаркированного
пикселя объекта запускает BFS (`BfsLabel`). Сложность по времени — O(rows × cols) на типичном изображении;
по памяти — O(rows × cols) для выходного буфера и очереди BFS.

## 4. Схема распараллеливания

Изображение делится на горизontальные полосы по числу потоков `num_threads = min(PPC_NUM_THREADS, rows)`.

| Этап                       | TBB-примитив                  | Диапазон                             | Синхронизация                  |
| -------------------------- | ----------------------------- | ------------------------------------ | ------------------------------ |
| Локальная маркировка полос | `task_arena` + `parallel_for` | `blocked_range<int>(0, num_threads)` | неявная в конце `parallel_for` |
| Смещение меток             | `parallel_for`                | по индексам полос                    | неявная                        |
| Union-find на границах     | последовательно               | границы полос                        | —                              |
| Первые позиции меток       | `task_arena` + `parallel_for` | потоки, затем корни                  | неявная                        |
| Финальный remap            | `task_arena` + `parallel_for` | `blocked_range<int>(0, rows)`        | неявная                        |

**Grainsize / partitioner:** используется `blocked_range` с параметрами по умолчанию (auto partitioner).
Для маркировки полос диапазон мал (`num_threads` ≤ 10), для remap — разбиение по строкам (2000 строк).

**Контроль конкуренции:** `tbb::task_arena arena(num_threads)` ограничивает число worker-потоков
значением `PPC_NUM_THREADS`.

**SEQ fallback:** при `num_threads <= 1` вызывается `RunSequentialMarking` без strip-пайплайна.

Ключевой фрагмент — маркировка полос:

```cpp
// File: tbb/src/ops_tbb.cpp
tbb::task_arena arena(num_threads);
arena.execute([&] {
  tbb::parallel_for(tbb::blocked_range<int>(0, num_threads),
                    [&](const tbb::blocked_range<int> &range) {
                      LabelStripsRange(input, output, cols, row_starts, labels_used,
                                       range.begin(), range.end());
                    });
});
```

Каждый subrange обрабатывает свои полосы независимо; метки пишутся сразу в `output` (без `local_planes`
на поток). Гонок на запись в `output` нет: полосы не пересекаются по строкам.

## 5. Детали реализации

Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.

Pipeline `RunImpl`:

1. `ParallelLabelStrips` — BFS внутри полосы (`BfsLabelInStrip`), очередь на `std::vector`.
2. `AddBaseOffsets` — добавление глобальных смещений меток.
3. `MergeBoundariesUnionFind` + `FlattenParentForest` — слияние компонент через границы полос.
4. `ComputeFirstPositionsParallel` + `BuildRemapFromFirstPositions` — перенумерация с сохранением порядка
   по первому вхождению метки.
5. `ApplyRemapInParallel` — remap только пикселей объектов.

`ValidationImpl`, `PreProcessingImpl`, `PostProcessingImpl` совпадают по контракту с SEQ.

## 6. Проверка корректности

Корректность TBB проверялась сравнением с SEQ на функциональных тестах (`tests/functional/main.cpp`):

| Тест                   | Описание                     |
| ---------------------- | ---------------------------- |
| `single_L_component`   | L-образная компонента 3×3    |
| `four_separate_pixels` | четыре изолированных объекта |
| `all_background`       | только фон                   |
| `all_objects`          | только объекты               |
| `two_horizontal_bars`  | две горизонтальные полосы    |

Результат: **5/5** тестов пройдено. Perf-тест дополнительно проверяет наличие меток > 0 на выходе.

## 7. Экспериментальная среда

| Параметр   | Значение                                                                  |
| ---------- | ------------------------------------------------------------------------- |
| CPU        | Apple M4, 10 ядер                                                         |
| ОС         | macOS 15.5                                                                |
| Компилятор | Apple Clang 17.0.0                                                        |
| Сборка     | Release, `USE_FUNC_TESTS=ON`, `USE_PERF_TESTS=ON`                         |
| Perf-вход  | 2000×2000, чётные строки — объект (`tests/performance/main.cpp`)          |
| Повторы    | 5 прогонов `RunImpl` на замер (инфраструктура курса), медиана по 3 сериям |
| Переменные | `PPC_NUM_THREADS` (экспортируется также как `OMP_NUM_THREADS`)            |

Команды:

```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
cmake --build build -j$(sysctl -n hw.ncpu)

./build/bin/ppc_func_tests --gtest_filter='*gaivoronskiy_m_marking_binary_components_tbb*'

PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/task_run_gaivoronskiy_m_marking_binary_components_tbb_enabled'
```

## 8. Результаты

Baseline SEQ — из `seq/report.md` (медиана 3 серий, perf-вход 2000×2000):

| Потоки   | SEQ, с   |
| -------- | -------- |
| 1        | 0.002901 |
| 2        | 0.002711 |
| 4        | 0.002683 |

Speedup = `T_seq / T_tbb` при том же числе потоков (см. корневой `report.md`).

### task_run

| Потоки   | TBB, с   | Speedup vs SEQ   | Efficiency   |
| -------- | -------- | ---------------- | ------------ |
| 1        | 0.003130 | 0.93             | 93%          |
| 2        | 0.001872 | 1.45             | 72%          |
| 4        | 0.001132 | **2.37**         | **59%**      |

### pipeline (4 потока)

| Backend   | Время, с   | Speedup vs SEQ   |
| --------- | ---------- | ---------------- |
| SEQ       | 0.007558   | 1.00             |
| TBB       | 0.003782   | 2.00             |

На 500×500 (ранние замеры) TBB не обгонял SEQ из‑за overhead планировщика. На 2000×2000 выигрыш
от параллелизации доминирует: ускорение ~2.4× при 4 потоках.

## 9. Выводы

TBB удобна для strip-маркировки: `parallel_for` по полосам и строкам, `task_arena` для контроля потоков.
Оптимизации (маркировка в `output`, `AddBaseOffsets`, параллельный remap) снижают память и число проходов.
На малых изображениях overhead TBB заметен; на 2000×2000 версия существенно быстрее SEQ и baseline TBB
без оптимизаций.
