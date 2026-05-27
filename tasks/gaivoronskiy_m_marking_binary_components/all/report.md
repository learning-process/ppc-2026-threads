# Маркировка компонент на бинарном изображении — ALL

- Student: Гайворонский Максим Витальевич, group 3823Б1ПР1
- Technology: ALL
- Variant: 1

## 1. Контекст

ALL-версия в рамках семестра «потоки» реализует оптимизированный strip-пайплайн маркировки компонент с
использованием `std::thread` (`TypeOfTask::kALL`). MPI-уровень в данной задаче не задействован: вычисления
выполняются в одном процессе с несколькими потоками. Последовательная версия (`seq/report.md`) — эталон
корректности.

## 2. Постановка задачи

Входной вектор: `[rows, cols, pixels...]`.

- `0` — пиксель объекта;
- `1` — фон.

Выходной вектор того же размера: фон = `0`, каждая 4-связная компонента — уникальная метка.

## 3. Базовый алгоритм

Последовательная BFS-маркировка (`seq/src/ops_seq.cpp`) — см. `seq/report.md`.

## 4. Межпроцессная схема

В текущей реализации **MPI не используется**: один процесс, `PPC_NUM_PROC = 1`. Конфигурация для замеров:
`1 process × N threads`, где `N = PPC_NUM_THREADS`.

При запуске под `mpirun` инфраструктура курса активирует func-тесты для `kALL`; локально без MPI func-тесты
пропускаются (см. раздел 6).

## 5. Внутрипроцессная схема

Реализация идентична оптимизированной STL-версии (`stl/src/ops_stl.cpp`): strip-декомпозиция + `std::thread`.

| Этап                  | Разбиение         | Потоки                          | Синхронизация   |
| --------------------- | ----------------- | ------------------------------- | --------------- |
| Маркировка полос      | по `tid`          | `num_threads - 1` worker + main | `join`          |
| AddBaseOffsets        | по `thread_idx`   | то же                           | `join`          |
| Union-find            | последовательно   | 1                               | —               |
| ComputeFirstPositions | по `tid` → строки | worker + main                   | `join`          |
| ApplyRemap            | по `tid` → строки | worker + main                   | `join`          |

**SEQ fallback** при `num_threads <= 1`.

Ключевой фрагмент:

```cpp
// File: all/src/ops_all.cpp
const int num_threads = NumThreadsForRows(rows);
if (num_threads <= 1) {
  RunSequentialMarking(input, output, rows, cols);
  return true;
}
ParallelLabelStrips(input, output, cols, num_threads, row_starts, labels_used);
AddBaseOffsets(output, cols, num_threads, row_starts, base);
// ... UF, remap ...
```

Гонок записи нет: потоки работают с непересекающимися диапазонами строк.

## 6. Детали реализации

Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.

Класс `GaivoronskiyMMarkingBinaryComponentsALL` возвращает `TypeOfTask::kALL`. Алгоритм `RunImpl` повторяет
оптимизированный STL-пайплайн:

1. Маркировка сразу в `output` (без `local_planes`).
2. `AddBaseOffsets` вместо копирования полос.
3. Union-find на границах полос.
4. Параллельный remap через `ComputeFirstPositionsParallel`.
5. BFS на `std::vector`.

## 7. Проверка корректности

| Способ                            | Результат                                        |
| --------------------------------- | ------------------------------------------------ |
| Func-тесты (`ppc_func_tests`)     | пропущены локально (требуют `mpirun` для `kALL`) |
| Perf-тест (`CheckTestOutputData`) | пройден — на выходе есть метки > 0               |
| Сравнение с SEQ/STL               | логика совпадает с оптимизированной STL-версией  |

Функциональные тесты (5 сценариев) зарегистрированы в `tests/functional/main.cpp` и выполняются под MPI-runner
курса.

## 8. Экспериментальная среда

| Параметр     | Значение                          |
| ------------ | --------------------------------- |
| CPU          | Apple M4, 10 ядер                 |
| ОС           | macOS 15.5                        |
| Компилятор   | Apple Clang 17.0.0, Release       |
| Perf-вход    | 2000×2000, чётные строки — объект |
| Конфигурация | 1 process × `PPC_NUM_THREADS`     |
| Повторы      | медиана по 3 сериям               |

Команды:

```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
cmake --build build -j$(sysctl -n hw.ncpu)

PPC_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/GaivoronskiyMMarkingPerfTests.RunPerfModes/task_run_gaivoronskiy_m_marking_binary_components_all_enabled'
```

## 9. Результаты

Baseline SEQ — из `seq/report.md` (медиана 3 серий, perf-вход 2000×2000):

| Потоки   | SEQ, с   |
| -------- | -------- |
| 1        | 0.002901 |
| 2        | 0.002711 |
| 4        | 0.002683 |

Speedup = `T_seq / T_all` при том же числе потоков (см. корневой `report.md`).

### task_run (1p × N threads)

| Потоки   | ALL, с   | Speedup vs SEQ   | Efficiency   |
| -------- | -------- | ---------------- | ------------ |
| 1        | 0.002640 | 1.10             | 110%         |
| 2        | 0.001846 | 1.47             | 73%          |
| 4        | 0.001021 | **2.63**         | **66%**      |

### pipeline (4 потока)

| Backend   | Время, с   | Speedup vs SEQ   |
| --------- | ---------- | ---------------- |
| SEQ       | 0.007558   | 1.00             |
| ALL       | 0.003799   | 1.99             |

На 2000×2000 ALL показывает ускорение **~2.6×** — сопоставимо со STL (~2.6×), поскольку вычислительное ядро
одинаково.

## 10. Выводы

ALL-версия в семестре потоков — оптимизированная thread-based реализация с типом `kALL` для инфраструктуры
курса. Без MPI-коммуникаций overhead минимален; результат близок к STL. Для полной проверки func-тестов
ALL требуется запуск под `mpirun` в CI.
