# Умножение плотных матриц блочным алгоритмом Кэннона

- Student: Синев Тимур Викторович, группа 3823Б1ПР5
- Variant: 1
- Local reports: [seq/report.md](seq/report.md), [omp/report.md](omp/report.md), [tbb/report.md](tbb/report.md),
  [stl/report.md](stl/report.md), [all/report.md](all/report.md)

## 1. Введение

В работе сравниваются **пять реализаций** умножения плотных матриц по
**блочному алгоритму Кэннона**: последовательная (SEQ) и параллельные — OpenMP
(OMP), oneTBB (TBB), `std::thread` (STL), гибрид MPI + OpenMP (ALL).

Задача удобна для курса по параллельному программированию: матрица делится на
блоки, на каждом шаге много **независимых** умножений блоков, между шагами —
только сдвиг блоков. Подробности по каждой технологии — в локальных отчётах;
здесь — общая постановка, методика, сводная таблица и выводы.

## 2. Единая постановка задачи

Даны две квадратные матрицы $A$ и $B$ размера $n \times n$ и размер блока
`b_size`. Нужно получить $C = A \cdot B$.

- **Вход:** `(b_size, A, B)` — см. `common/include/common.hpp`.
- **Выход:** матрица $C$ того же размера.
- **Ограничения:** $n$ делится на `b_size` без остатка; матрицы квадратные.
- **Критерий корректности:** совпадение с эталоном из функциональных тестов,
  погрешность по элементам не больше $10^{-10}$.

Алгоритм Кэннона: сетка блоков $grid\_sz = n / b\_size$, на каждом из
$grid\_sz$ шагов — умножение и накопление в блоках результата, затем
циклический сдвиг блоков $A$ и $B$.

## 3. Единая методика эксперимента

**Окружение:**

- **CPU:** AMD RYZEN 5 2600.
- **RAM:** 16 GB.
- **ОС:** Windows 11.
- **Build type:** Release.
- **Compiler:** Microsoft Visual C++ 19.44 (MSVC 2022).

**Методика:**

- Размер в perf-тесте: матрицы **$512 \times 512$**, блок **$32 \times 32$**
  (`tests/performance/main.cpp`).
- Для OMP, TBB, STL число потоков: `PPC_NUM_THREADS` = 1, 2, 4, 8
  (также `OMP_NUM_THREADS`).
- Для ALL — запуск через `mpiexec`; всего **$p$** работников =
  процессы × потоки (см. [all/report.md](all/report.md)).
- Замеры в режиме **`task_run`** (основное вычисление в `RunImpl`).
- Базовое время: **$T_{seq} = 0.0499$ с** — SEQ при **одном потоке**.
- Ускорение: $S = T_{seq} / T$.
- Эффективность: $E = S / p \cdot 100\%$, где $p$ — число потоков
  (или работников для ALL).

Дата замеров: **18.05.2026**, по одному прогону `ppc_perf_tests`
на каждую конфигурацию.

## 4. Сводка корректности

Все реализации прошли **8 функциональных тестов** (размеры от $2 \times 2$ до
$9 \times 9$, разные `b_size`):

| Backend | Результат |
| ------- | --------- |
| SEQ, OMP, TBB, STL | 8/8 PASSED |
| ALL | 8/8 PASSED под `mpiexec -n 2` |

Результаты совпадают с эталонными матрицами из `tests/functional/main.cpp`,
значит параллельные версии считают то же, что и SEQ.

## 5. Агрегированные результаты

В таблице — **лучшие** показатели по каждому backend-у в наших замерах
(для параллельных версий — при максимальном $p = 8$). Полные таблицы по всем
$p$ — в локальных отчётах.

| Backend | Конфигурация | Время (с) | Ускорение ($S$) | Эффективность ($E$) | Примечание |
| ------- | ------------ | --------- | --------------- | ------------------- | ---------- |
| **SEQ** | 1 поток | 0.0499 | 1.00 | 100% | Baseline |
| **OMP** | 8 потоков | 0.0238 | 2.10 | 26.2% | Стабильное ускорение |
| **TBB** | 8 потоков | 0.0160 | **3.12** | 39.0% | **Самая быстрая версия** |
| **STL** | 8 потоков | 0.0533 | 0.94 | 11.7% | Медленнее baseline |
| **ALL** | 2 процесса × 4 потока | 0.0253 | 1.98 | 24.7% | MPI + OpenMP |

### Динамика по числу потоков (режим `task_run`)

Кратко по локальным отчётам:

- **TBB** — время падает при росте $p$ (0.0585 → 0.0160 с);
  лучший рост ускорения.
- **OMP** — тоже ускоряется (до ~2.1× на 8 потоках), но чуть слабее TBB.
- **STL** — на 1 потоке сильно медленнее SEQ; на 4 потоках почти как SEQ;
  на 8 снова хуже.
- **SEQ** — при $p > 1$ время почти не меняется (код не параллельный,
  разброс — шум замеров).
- **ALL** — на 1 процессе медленнее SEQ; с ростом $p$ становится быстрее.

## 6. Интерпретация различий

1. **SEQ (эталон).** Задаёт базовое время $T_{seq} = 0.0499$ с.
   Основная работа — умножение блоков в тройных циклах и $grid\_sz$ шагов
   Кэннона. Эту версию используем для проверки правильности остальных.

2. **OpenMP.** Циклы по блокам `(i, j)` распараллены через
   `#pragma omp parallel for`. На 8 потоках ускорение **~2.1×**.
   Ограничение — барьер в конце каждой параллельной области
   и шаги сдвига блоков, которые выполняются после синхронизации.

3. **oneTBB (лучший результат).** Сетка блоков обрабатывается через
   `tbb::parallel_for` и `blocked_range2d`. На 8 потоках **~3.1×** к SEQ.
   Библиотека сама делит работу на части и перераспределяет нагрузку
   между потоками — для равномерных блоков это выгодно.

4. **STL (`std::thread`).** Потоки создаются вручную в `ParallelFor`,
   в конце — `join`. На 1 потоке версия **вдвое медленнее** SEQ
   из‑за лишних затрат на потоки. При $grid\_sz = 16$ мало крупных частей
   работы — ручные потоки проигрывают OMP и TBB.

5. **ALL (MPI + OpenMP).** Строки блоков делятся между MPI-процессами,
   внутри процесса — OpenMP. Нужен `mpiexec`. На одном процессе медленнее
   SEQ (расходы MPI). При **2 × 4** работниках ускорение **~2×**,
   но чистый OMP/TBB на одной машине всё же быстрее из‑за обмена данными
   (`Bcast`, `Allgatherv`).

## 7. Репродуцируемость

```text
# Сборка
cmake -S . -B build
cmake --build build --config Release --target ppc_func_tests ppc_perf_tests

# Функциональные тесты
build\bin\ppc_func_tests.exe --gtest_filter=*timur_a_cannon*

# Perf: SEQ (пример для 1 потока)
set PPC_NUM_THREADS=1
set OMP_NUM_THREADS=1
set PPC_NUM_PROC=1
build\bin\ppc_perf_tests.exe --gtest_filter=*timur_a_cannon_seq*

# Perf: OMP / TBB / STL (пример для 8 потоков)
set PPC_NUM_THREADS=8
set OMP_NUM_THREADS=8
build\bin\ppc_perf_tests.exe --gtest_filter=*timur_a_cannon_omp*

# Perf: ALL (2 процесса, 4 потока в каждом — всего 8 работников)
mpiexec -env PPC_NUM_THREADS 4 -env OMP_NUM_THREADS 4 \
  -env PPC_NUM_PROC 2 -n 2 \
  build\bin\ppc_perf_tests.exe --gtest_filter=*timur_a_cannon_all*

Альтернатива из курса: `scripts/run_tests.py --running-type=threads --counts 1 2 4 8`
и `--running-type=processes` для ALL.

## 8. Заключение

Для умножения матриц **$512 \times 512$** алгоритмом Кэннона на тестовой машине **самой быстрой**
оказалась реализация на **oneTBB** (ускорение **~3.1×** к SEQ на 8 потоках).
**OpenMP** — хороший практичный вариант (**~2.1×**). 
**STL** без пула потоков на этой задаче **не даёт выигрыша**. 
**Гибрид MPI + OpenMP** имеет смысл при нескольких процессах и больших матрицах; на малых размерах overhead MPI заметен.

Для защиты работы детальные таблицы, фрагменты кода и комментарии по каждой технологии см. в
[seq/report.md](seq/report.md), [omp/report.md](omp/report.md), [tbb/report.md](tbb/report.md),
[stl/report.md](stl/report.md), [all/report.md](all/report.md).

## 9. Источники

1. Материалы курса «Параллельное программирование», репозиторий
   [ppc-2026-threads](https://github.com/learning-process/ppc-2026-threads).
2. [OpenMP](https://www.openmp.org/) — директивы `parallel for`, `collapse`.
3. [oneTBB](https://github.com/uxlfoundation/oneTBB) — `parallel_for`, `blocked_range2d`.
4. [MPI Forum](https://www.mpi-forum.org/) — `MPI_Bcast`, `MPI_Allgatherv`.
5. [cppreference.com](https://en.cppreference.com/) — `std::thread`, работа с потоками.

## 10. Приложение

Короткие листинги, дополнительные графики, поясняющие диаграммы.

### Структура задачи

Каталог задачи построен по минимальному каркасу курса

```text
tasks/timur_a_cannon/
  report.md                      # обязательный корневой сводный отчёт
  info.json                      # сведения о студенте
  settings.json                  # включённые технологии
  common/
    include/common.hpp           # InType, OutType, TestType, BaseTask
  seq/
    include/ops_seq.hpp
    src/ops_seq.cpp
    report.md                    # локальный отчёт по SEQ
  omp/
    include/ops_omp.hpp
    src/ops_omp.cpp
    report.md                    # локальный отчёт по OMP
  tbb/
    include/ops_tbb.hpp
    src/ops_tbb.cpp
    report.md                    # локальный отчёт по TBB
  stl/
    include/ops_stl.hpp
    src/ops_stl.cpp
    report.md                    # локальный отчёт по std::thread
  all/
    include/ops_all.hpp
    src/ops_all.cpp
    report.md                    # локальный отчёт по гибридной версии
  tests/
    functional/main.cpp
    performance/main.cpp
  data/                          # опционально (в задаче не используется)
  img/                           # опционально 
```

В `common/include/common.hpp` заданы типы задачи `InType`, `OutType`,
`TestType`, `BaseTask`. В каждом каталоге `seq/`, `omp/`, `tbb/`, `stl/`,
`all/` — класс-наследник `BaseTask` со своим `TypeOfTask` и методами
`ValidationImpl`, `PreProcessingImpl`, `RunImpl`, `PostProcessingImpl`.

`tests/functional/main.cpp` — один набор из **восьми** тестовых случаев
(`a`–`h`, размеры от $2 \times 2$ до $9 \times 9$) для всех backend-ов;
эталон — заранее вычисленная матрица произведения.

`tests/performance/main.cpp` — общий каркас курса (`BaseRunPerfTests`,
`MakeAllPerfTasks`), режимы `task_run` и `pipeline`, матрицы $512 \times 512$,
`b_size = 32`.

Каталог `all/` — гибридная версия: `MPI_Comm_rank`, `MPI_Bcast` входных
матриц, распределение строк блоков между процессами, `MPI_Allgatherv`
результата, внутри процесса — OpenMP в `ComputeLocalResult`
(подробнее — [all/report.md](all/report.md)).

### Короткий листинг: типы входа и выхода

```cpp
// File: common/include/common.hpp
using InType = std::tuple<int, std::vector<std::vector<double>>,
                          std::vector<std::vector<double>>>;
using OutType = std::vector<std::vector<double>>;
using TestType = std::tuple<std::string, int,
                            std::vector<std::vector<double>>,
                            std::vector<std::vector<double>>,
                            std::vector<std::vector<double>>>;
using BaseTask = ppc::task::Task<InType, OutType>;
```

`InType`: размер блока `b_size`, матрицы `A` и `B`.  
`OutType`: матрица `C = A·B`.
