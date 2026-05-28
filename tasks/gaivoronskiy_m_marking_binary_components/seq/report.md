# Маркировка компонент на бинарном изображении — SEQ

- Student: Гайворонский Максим Витальевич, group 3823Б1ПР1
- Technology: SEQ
- Variant: 1

## 1. Контекст

Последовательная реализация выполняет маркировку 4-связных компонент объектов на бинарном изображении.
Она служит **эталоном корректности** для параллельных backend-ов (OMP, STL, TBB, ALL) и **baseline**
для расчёта speedup и efficiency в сводном отчёте (`../report.md`).

Код SEQ одинаков во всех ветках задачи и не содержит strip-декомпозиции или union-find — только
построчный обход и BFS по каждой новой компоненте.

## 2. Постановка задачи

Входной вектор: `[rows, cols, pixels...]`.

| Значение   | Смысл           |
| ---------- | --------------- |
| `0`        | пиксель объекта |
| `1`        | фон             |

Выходной вектор того же размера: фон = `0`, каждая 4-связная компонента — уникальная положительная метка
(`1`, `2`, …). Связность — только по четырём соседям (верх, низ, лево, право).

Ограничения проверяются в `ValidationImpl`: `rows > 0`, `cols > 0`, размер вектора равен `rows * cols + 2`.

## 3. Алгоритм

1. **PreProcessingImpl** — выделить выходной буфер, записать `rows` и `cols` в первые два элемента.
2. **RunImpl** — обход изображения построчно слева направо, сверху вниз.
3. Для каждого непромаркированного пикселя объекта (`input[idx] == 0 && output[idx] == 0`):
   - увеличить счётчик меток;
   - запустить **BFS** (`BfsLabel`) от этой точки.
4. BFS помещает стартовую координату в очередь, присваивает метку и распространяет её на все
   4-соседние пиксели объекта, ещё не получившие метку.
5. **PostProcessingImpl** — без дополнительных действий.

**Сложность по времени:** O(rows × cols) — каждый пиксель объекта посещается BFS не более одного раза.

**Сложность по памяти:** O(rows × cols) для выходного буфера плюс O(размер компоненты) для очереди BFS
(в худшем случае O(rows × cols)).

## 4. Детали реализации

Файлы: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`.

| Метод                | Назначение                                                    |
| -------------------- | ------------------------------------------------------------- |
| `ValidationImpl`     | проверка размеров и согласованности длины вектора             |
| `PreProcessingImpl`  | инициализация `GetOutput()` нулями, копирование `rows`/`cols` |
| `RunImpl`            | двойной цикл + вызовы `BfsLabel`                              |
| `PostProcessingImpl` | всегда `true`                                                 |

Ключевой фрагмент — запуск BFS для новой компоненты:

```cpp
// File: seq/src/ops_seq.cpp
for (int i = 0; i < rows; i++) {
  for (int j = 0; j < cols; j++) {
    int idx = (i * cols) + j + 2;
    if (input[idx] == 0 && output[idx] == 0) {
      label++;
      BfsLabel(input, output, rows, cols, i, j, label);
    }
  }
}
```

Внутри `BfsLabel` используется `std::queue<std::pair<int,int>>` и массивы смещений `kDx`/`kDy` для
4-соседства. Метки записываются сразу в `output[idx]`; повторное посещение отсекается условием
`output[nidx] == 0`.

Переменная окружения `PPC_NUM_THREADS` на алгоритм **не влияет** — вычисление полностью однопоточное.

## 5. Проверка корректности

Общий набор функциональных тестов (`tests/functional/main.cpp`) прогоняет все backend-ы на одних данных.

| Тест                   | Описание                     |
| ---------------------- | ---------------------------- |
| `single_L_component`   | L-образная компонента 3×3    |
| `four_separate_pixels` | четыре изолированных объекта |
| `all_background`       | только фон                   |
| `all_objects`          | только объекты               |
| `two_horizontal_bars`  | две горизонтальные полосы    |

Результат для SEQ: **5/5** функциональных тестов пройдено.

Параллельные реализации сверяются с SEQ побайтово на тех же входах.

## 6. Экспериментальная среда

| Параметр    | Значение                                                             |
| ----------- | -------------------------------------------------------------------- |
| CPU         | Apple M4, 10 ядер                                                    |
| ОС          | macOS 15.5                                                           |
| Компилятор  | Apple Clang 17.0.0                                                   |
| Сборка      | Release, `USE_FUNC_TESTS=ON`, `USE_PERF_TESTS=ON`                    |
| Perf-вход   | **2000×2000**, чётные строки — объект (`tests/performance/main.cpp`) |
| Повторы     | медиана по 3 сериям                                                  |
| Режимы perf | `task_run` и `pipeline`                                              |
| Переменные  | `PPC_NUM_THREADS` (для SEQ не меняет алгоритм)                       |

Команды:

```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
cmake --build build -j$(sysctl -n hw.ncpu)

./build/bin/ppc_func_tests --gtest_filter='*gaivoronskiy_m_marking_binary_components_seq*'

PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/task_run_gaivoronskiy_m_marking_binary_components_seq_enabled'
```

## 7. Результаты

SEQ — baseline: speedup = 1.00, efficiency не определена (один рабочий поток).

### task_run (медиана 3 серий, 2000×2000)

| Потоки   | Время, с   | Speedup   | Efficiency   |
| -------- | ---------- | --------- | ------------ |
| 1        | 0.002901   | 1.00      | N/A          |
| 2        | 0.002711   | 1.00      | N/A          |
| 4        | 0.002683   | 1.00      | N/A          |

Значения при 2 и 4 «потоках» совпадают с 1-поточным режимом с точностью до шума измерений —
алгоритм не использует параллелизм.

### pipeline (4 потока, 2000×2000)

| Режим    | Время, с   |
| -------- | ---------- |
| pipeline | 0.007558   |

Режим `pipeline` включает полный цикл `ValidationImpl` → `PreProcessingImpl` → `RunImpl` →
`PostProcessingImpl`; для SEQ основное время даёт `RunImpl` (~2.7 ms), остальные этапы добавляют
накладные расходы.

На размере **500×500** SEQ занимает ~0.15 ms и часто быстрее параллельных версий из‑за overhead
strip-пайплайна. На **2000×2000** (~4M пикселей) baseline составляет ~2.7 ms — здесь параллельные
backend-ы получают устойчивое ускорение (см. `omp/report.md`, `stl/report.md`, `tbb/report.md`,
`all/report.md` и корневой `report.md`).

## 8. Выводы

Последовательный BFS прост в реализации, детерминирован и полностью покрывает функциональные сценарии.
Он задаёт эталонные метки компонент и baseline для сравнения технологий.

Ограничения SEQ очевидны: один поток, последовательный обход, очередь BFS с произвольными обращениями
к памяти. Параллельные реализации компенсируют это strip-декомпозицией и union-find на границах полос,
но выигрывают только когда стоимость вычисления превышает накладные расходы синхронизации.

## 9. Репродуцируемость

Фильтры perf-тестов:

```bash
# task_run
PPC_NUM_THREADS=1 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/task_run_gaivoronskiy_m_marking_binary_components_seq_enabled'

# pipeline
PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/pipeline_gaivoronskiy_m_marking_binary_components_seq_enabled'
```

Сводное сравнение backend-ов — в `../report.md`.
