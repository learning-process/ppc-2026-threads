# Обработка контуров бинарных компонент — Сводный отчёт

- **Student:** Перяшкин Василий Андреевич, группа 3823Б1ПР4
- **Variant:** 30
- **Локальные отчёты:** [seq/report.md](seq/report.md), [omp/report.md](omp/report.md), [tbb/report.md](tbb/report.md), [stl/report.md](stl/report.md), [all/report.md](all/report.md)

## 1. Введение

Задача — обработка бинарного изображения: нахождение всех 4-связных компонент и вычисление выпуклой оболочки каждой из них. Алгоритм естественно делится на два этапа: последовательный BFS и embarrassingly parallel вычисление оболочек. Это делает задачу подходящей для сравнения нескольких моделей параллелизма: все они применяют параллелизм ровно к той же фазе, но разными средствами.

## 2. Единая постановка задачи

**Вход:** `BinaryImage{width, height, data}` — бинарный растр (0/1).  
**Выход:** `vector<vector<Point>>` — по одной выпуклой оболочке на компоненту.  
**Критерий корректности:** результат идентичен SEQ-эталону побайтно (порядок компонент и вершин детерминирован).  
**Граничные случаи:** пустое изображение, одиночный пиксель, линия, прямоугольник, касающийся границы.

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

**Тестовый набор для замеров:** бинарное изображение 512×512 с шахматным паттерном (шаг 2). Генерация:
```cpp
// tests/performance/main.cpp — MakePattern(512, 512, 2)
```

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

**Ограничения применимости STL:** реализация использует `std::ranges::for_each` без политики выполнения — фактически последовательна и не масштабируется с числом потоков.

## 5. Агрегированные результаты

### Режим task

| Backend | Workers     | Медианное время (task) | Speedup (vs SEQ) | Efficiency |
|---------|-------------|------------------------|------------------|------------|
| SEQ     | 1           | 1.370 мс               | 1.00×            | 100%       |
| OMP     | 1           | 1.316 мс               | 1.04×            | 104%       |
| OMP     | 2           | 1.285 мс               | 1.07×            | 53%        |
| OMP     | 4           | 1.320 мс               | 1.04×            | 26%        |
| TBB     | 1           | 1.402 мс               | 0.98×            | 98%        |
| TBB     | 2           | 1.036 мс               | 1.32×            | 66%        |
| TBB     | 4           | 1.032 мс               | 1.33×            | 33%        |
| STL     | 1 (факт.)   | 1.364 мс               | 1.00×            | —          |
| ALL     | 2×1         | 1.379 мс               | 0.99×            | 50%        |
| ALL     | 4×1         | 2.565 мс               | 0.53×            | 13%        |
| ALL     | 2×2         | 1.393 мс               | 0.98×            | 25%        |

### Режим pipeline

| Backend | Workers     | Медианное время (pipeline) | Speedup (vs SEQ) | Efficiency |
|---------|-------------|----------------------------|------------------|------------|
| SEQ     | 1           | 1.452 мс                   | 1.00×            | 100%       |
| OMP     | 4           | 1.252 мс                   | 1.16×            | 29%        |
| TBB     | 4           | 0.910 мс                   | 1.60×            | 40%        |
| STL     | 1 (факт.)   | 1.359 мс                   | 1.07×            | —          |
| ALL     | 2×2         | 1.472 мс                   | 0.99×            | 25%        |

## 6. Интерпретация различий

**SEQ** — чистый baseline. Самый дорогой фрагмент — BFS (~60–70% времени), который не распараллелен ни в одной версии.

**OMP** — минимальные изменения кода, ускорение ограничено законом Амдала. При доле параллельного участка ~35% теоретический предел для 4 потоков составляет ≈1.4–1.5×. `schedule(static)` подходит для равномерной нагрузки (шахматный паттерн).

**TBB** — семантически идентичен OMP-версии: тот же диапазон задач, те же граничные условия. `auto_partitioner` адаптивно настраивает grain size — преимущество проявится при неравномерных компонентах. На шахматном паттерне результаты OMP и TBB ожидаемо близки, что объясняется схожей grain-структурой работы и одинаковым узким местом (BFS).

**STL** — реализован через `std::ranges::for_each` без политики выполнения `std::execution::par`, поэтому фактически последователен. Время совпадает с SEQ. Для настоящего параллелизма потребовалась бы замена на `std::for_each(std::execution::par_unseq, ...)`.

**ALL** — при малых размерах задачи проигрывает SEQ из-за overhead на сериализацию компонент и MPI-коммуникации. Преимущество реализуется только при большом числе компонент (сотни тысяч), когда стоимость MPI-вызовов амортизируется. `MPI_Bcast` финального результата — потенциальная точка насыщения при большом числе рангов.

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

Лучшая версия для одной машины с разделяемой памятью — **TBB** при числе потоков 4 (pipeline: 0.910 мс, speedup 1.60×). OMP даёт скромное ускорение 1.04–1.07× из-за доминирования последовательного BFS. STL-версия не масштабируется в текущей реализации (фактически SEQ). Гибридная ALL-версия **проигрывает** SEQ на малых задачах: overhead MPI-сериализации и коммуникаций превышает выигрыш от распределения вычислений — 0.53–0.99× против SEQ на 512×512. ALL оправдан только при больших изображениях с сотнями тысяч компонент.

Ограничения сравнения: все замеры выполнены в контейнере, что может вносить шум из-за ограниченных аппаратных ресурсов. Числа нельзя механически переносить на другое железо.

Что можно улучшить:
- Реализовать параллельный BFS (например, spanning-tree подход с атомарными операциями на `vis`).
- В STL заменить `std::ranges::for_each` на `std::for_each(std::execution::par_unseq, ...)`.
- В ALL убрать финальный `MPI_Bcast` для версий, где ранг 0 возвращает результат (оптимизация для однопользовательского сценария).

## 9. Источники

- Документация курса (tasks/example_threads)
- Спецификация OpenMP 5.2: [openmp.org/specifications](https://www.openmp.org/specifications/)
- Документация oneTBB от UXL Foundation: [uxlfoundation.github.io/oneTBB](https://uxlfoundation.github.io/oneTBB/)
- MPI Forum, MPI-4.0 Standard: [mpi-forum.org/docs](https://www.mpi-forum.org/docs/)
- cppreference.com: `std::ranges::for_each`, `std::execution::par`

## 10. Приложение

### Схема алгоритма SEQ

```
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

```
Rank 0: BFS → flatten → MPI_Scatterv →→→→→→→→→→→→→→→→→→→
                                        Rank i: unflatten → OMP hulls → flatten
Rank 0: MPI_Gatherv ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
        unflatten → MPI_Bcast →→→→→→→→→→→→→→→→→→→→→→→→→→→
                                        Rank i: unflatten → local_out_
```
