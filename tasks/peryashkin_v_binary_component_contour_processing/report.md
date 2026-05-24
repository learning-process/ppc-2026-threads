# Обработка контуров бинарных компонент — Сводный отчёт

- **Student:** Перяшкин Василий Андреевич, группа 3823Б1ПР4
- **Variant:** 30
- **Локальные отчёты:** [seq](seq/report.md), [omp](omp/report.md), [tbb](tbb/report.md),
  [stl](stl/report.md), [all](all/report.md)

## 1. Введение

Задача — обработка бинарного изображения: нахождение всех 4-связных компонент и вычисление
выпуклой оболочки каждой из них. Алгоритм делится на два этапа: последовательный BFS и
embarrassingly parallel вычисление оболочек. Все backend'ы применяют параллелизм к той же
фазе, но разными средствами.

## 2. Единая постановка задачи

**Вход:** `BinaryImage{width, height, data}` — бинарный растр (0/1).
**Выход:** `vector<vector<Point>>` — по одной выпуклой оболочке на компоненту.
**Критерий корректности:** результат идентичен SEQ-эталону побайтно.
**Граничные случаи:** пустое изображение, одиночный пиксель, линия, прямоугольник у границы.

## 3. Единая методика эксперимента

| Параметр         | Значение                                      |
|------------------|-----------------------------------------------|
| CPU              | Apple M2 (4P + 4E cores, ARM64)               |
| Ядра / потоки    | 8 / 8                                         |
| RAM              | 8 GB                                          |
| OS               | macOS 13″ / Docker (Linux container)          |
| Компилятор       | GCC 13.3.0, `-std=c++20`, `-O3`               |
| MPI              | Open MPI 4.1.6                                |
| TBB              | oneTBB 2022.0.0                               |
| Build type       | Release (`-O3 -DNDEBUG`)                      |

**Тестовый набор для замеров:** бинарное изображение 512×512 с шахматным паттерном (шаг 2).
Генерация: `MakePattern(512, 512, 2)` в `tests/performance/main.cpp`.

**Определения метрик (единые для всех таблиц):**

- `speedup = T_seq / T_backend`
- `efficiency = speedup / workers × 100%`
- `workers` = число потоков (OMP/TBB/STL) или `ranks × threads_per_rank` (ALL)
- Медиана по 5 повторным запускам

**Команды сборки:**

```bash
git submodule update --init --recursive --depth=1
cmake -S . -B build -DUSE_PERF_TESTS=ON -DUSE_FUNC_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

**Команды запуска функциональных тестов:**

```bash
PPC_NUM_THREADS=4 scripts/run_tests.py --running-type=threads --counts 1 2 4 \
  --gtest_filter="*peryashkin_v*"
PPC_NUM_PROC=2 PPC_NUM_THREADS=2 scripts/run_tests.py --running-type=processes \
  --gtest_filter="*peryashkin_v*all*"
```

**Команды замеров производительности:**

```bash
PPC_NUM_THREADS=4 scripts/run_tests.py --running-type=performance \
  --gtest_filter="*peryashkin_v*"
```

## 4. Сводка корректности

Все пять backend'ов дают результат, идентичный SEQ, на всех 7 функциональных тестах:

| Backend | empty | point | line | square | two_components | hole | touch_border |
|---------|-------|-------|------|--------|----------------|------|--------------|
| SEQ     | ✓     | ✓     | ✓    | ✓      | ✓              | ✓    | ✓            |
| OMP     | ✓     | ✓     | ✓    | ✓      | ✓              | ✓    | ✓            |
| TBB     | ✓     | ✓     | ✓    | ✓      | ✓              | ✓    | ✓            |
| STL     | ✓     | ✓     | ✓    | ✓      | ✓              | ✓    | ✓            |
| ALL     | ✓     | ✓     | ✓    | ✓      | ✓              | ✓    | ✓            |

**Ограничения STL:** реализация использует `std::ranges::for_each` без политики выполнения —
фактически последовательна и не масштабируется с числом потоков.

## 5. Агрегированные результаты

### Режим task

| Backend | Workers   | Медианное время (task) | Speedup (vs SEQ) | Efficiency |
|---------|-----------|------------------------|------------------|------------|
| SEQ     | 1         | 1.370 мс               | 1.00×            | 100%       |
| OMP     | 1         | 1.316 мс               | 1.04×            | 104%       |
| OMP     | 2         | 1.285 мс               | 1.07×            | 53%        |
| OMP     | 4         | 1.320 мс               | 1.04×            | 26%        |
| TBB     | 1         | 1.402 мс               | 0.98×            | 98%        |
| TBB     | 2         | 1.036 мс               | 1.32×            | 66%        |
| TBB     | 4         | 1.032 мс               | 1.33×            | 33%        |
| STL     | 1 (факт.) | 1.364 мс               | 1.00×            | —          |
| ALL     | 2×1       | 1.379 мс               | 0.99×            | 50%        |
| ALL     | 4×1       | 2.565 мс               | 0.53×            | 13%        |
| ALL     | 2×2       | 1.393 мс               | 0.98×            | 25%        |

### Режим pipeline

| Backend | Workers   | Медианное время (pipeline) | Speedup (vs SEQ) | Efficiency |
|---------|-----------|----------------------------|------------------|------------|
| SEQ     | 1         | 1.452 мс                   | 1.00×            | 100%       |
| OMP     | 4         | 1.252 мс                   | 1.16×            | 29%        |
| TBB     | 4         | 0.910 мс                   | 1.60×            | 40%        |
| STL     | 1 (факт.) | 1.359 мс                   | 1.07×            | —          |
| ALL     | 2×2       | 1.472 мс                   | 0.99×            | 25%        |

## 6. Интерпретация различий

**SEQ** — чистый baseline. Самый дорогой фрагмент — BFS (~60–70% времени),
который не распараллелен ни в одной версии.

**OMP** — минимальные изменения кода, ускорение ограничено законом Амдала.
При доле параллельного участка ~35% теоретический предел для 4 потоков ≈ 1.4–1.5×.
`schedule(static)` подходит для равномерной нагрузки (шахматный паттерн).

**TBB** — семантически идентичен OMP-версии: тот же диапазон задач, те же граничные условия.
`auto_partitioner` адаптивно настраивает grain size — преимущество при неравномерных компонентах.
На шахматном паттерне результаты OMP и TBB близки из-за одинакового узкого места (BFS).

**STL** — реализован через `std::ranges::for_each` без политики `std::execution::par`,
поэтому фактически последователен. Для настоящего параллелизма нужна замена на
`std::for_each(std::execution::par_unseq, ...)`.

**ALL** — при малых размерах задачи проигрывает SEQ из-за overhead на сериализацию и
MPI-коммуникации. Преимущество реализуется только при большом числе компонент (сотни тысяч),
когда стоимость MPI-вызовов амортизируется.

## 7. Воспроизводимость

```bash
# Сборка
cmake -S . -B build -DUSE_PERF_TESTS=ON -DUSE_FUNC_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# Функциональные тесты (все backend'ы)
PPC_NUM_THREADS=4 scripts/run_tests.py --running-type=threads --counts 1 2 4 \
  --gtest_filter="*peryashkin_v*"
PPC_NUM_PROC=4 PPC_NUM_THREADS=1 scripts/run_tests.py --running-type=processes \
  --gtest_filter="*peryashkin_v*all*"

# Замеры производительности
PPC_NUM_THREADS=4 scripts/run_tests.py --running-type=performance \
  --gtest_filter="*peryashkin_v*"
PPC_NUM_PROC=4 scripts/run_tests.py --running-type=performance \
  --gtest_filter="*peryashkin_v*all*"
```

## 8. Заключение

Лучшая версия для одной машины с разделяемой памятью — **TBB** при 4 потоках
(pipeline: 0.910 мс, speedup 1.60×). OMP даёт скромное ускорение 1.04–1.07× из-за
доминирования последовательного BFS. STL-версия не масштабируется в текущей реализации.
Гибридная ALL-версия **проигрывает** SEQ на малых задачах: overhead MPI-сериализации
превышает выигрыш от распределения вычислений (0.53–0.99× на 512×512).
ALL оправдан только при больших изображениях с сотнями тысяч компонент.

Ограничения: замеры выполнены в контейнере — числа нельзя механически переносить на другое железо.

Что можно улучшить:

- Реализовать параллельный BFS (spanning-tree подход с атомарными операциями на `vis`).
- В STL заменить `std::ranges::for_each` на `std::for_each(std::execution::par_unseq, ...)`.
- В ALL убрать финальный `MPI_Bcast` там, где ранг 0 достаточен для возврата результата.

## 9. Источники

- Документация курса (`tasks/example_threads`)
- Спецификация OpenMP 5.2: [openmp.org/specifications](https://www.openmp.org/specifications/)
- Документация oneTBB: [uxlfoundation.github.io/oneTBB](https://uxlfoundation.github.io/oneTBB/)
- MPI Forum, MPI-4.0 Standard: [mpi-forum.org/docs](https://www.mpi-forum.org/docs/)
- cppreference.com: `std::ranges::for_each`, `std::execution::par`

## 10. Приложение

### Схема алгоритма SEQ

```text
BinaryImage (W×H)
       │
       ▼
ExtractComponents4   ← BFS, 4-связность, последовательно, O(W×H)
       │
       ▼
[comp_0, comp_1, ..., comp_K]
       │  ← параллелизуется в OMP/TBB/ALL
       ▼
ConvexHullMonotonicChain × K   ← O(M_k log M_k) на компоненту
       │
       ▼
[hull_0, hull_1, ..., hull_K]
```

### Схема MPI-обмена в ALL

```text
Rank 0: BFS → flatten → MPI_Scatterv →→→→→→→→→→→→→→→→→→→
                                        Rank i: unflatten → OMP hulls → flatten
Rank 0: MPI_Gatherv ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
        unflatten → MPI_Bcast →→→→→→→→→→→→→→→→→→→→→→→→→→→
                                        Rank i: unflatten → local_out_
```
