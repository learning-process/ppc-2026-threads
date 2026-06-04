# Поразрядная сортировка double с четно-нечетным слиянием Бэтчера

- Студент: Гусев Дмитрий Алексеевич
- Группа: 3823Б1ФИ1
- Вариант: 20
- Технологии: SEQ, OMP, TBB, STL, ALL
- Локальные отчеты: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

---

## 1. Введение

В работе исследуется поразрядная сортировка вещественных чисел типа `double` с
четно-нечетным слиянием Бэтчера. Задача реализована в пяти вариантах: последовательном
SEQ, OpenMP, oneTBB, `std::thread` и гибридном ALL.

Алгоритм удобен для сравнения разных моделей параллелизма, потому что имеет естественную
блочную структуру. Сначала независимые части массива можно сортировать отдельно, затем
отсортированные блоки объединяются раундами слияния. При этом корректность для `double`
зависит от правильного преобразования битового представления в сортируемый ключ.

---

## 2. Единая постановка задачи

Вход: `std::vector<double>`.

Выход: вектор той же длины, отсортированный по неубыванию.

Общие критерии корректности:

1. Результат совпадает с `std::ranges::sort`.
2. Сохраняется кратность элементов.
3. Корректно обрабатываются пустые и малые входы.
4. Поддерживаются отрицательные значения, положительные значения, экстремальные числа и
   знаковые нули.
5. Параллельные реализации не имеют гонок данных на выходном массиве.

---

## 3. Единая методика эксперимента

### Окружение

- Процессор: Intel(R) Core(TM) Ultra 9 185H
- Оперативная память: 31.4 GB
- ОС: Windows 11
- Компилятор: MSVC 19.44.35220
- CMake: 4.2.0-rc2
- Тип сборки: Release

### Переменные и метрики

- `PPC_NUM_THREADS` задает число потоков для OMP, TBB, STL и внутрипроцессной части ALL.
- `PPC_NUM_PROC` и `mpiexec -n` задают число процессов для ALL.
- `time` - измеренное время `Run` или `task_run`.
- `speedup = baseline_time / current_time`.
- `efficiency = speedup / workers`.

Для SEQ `workers = 1`. Для OMP и TBB в локальных отчетах дополнительно приведена
таблица перебора `workers`. Для STL стандартный тест производительности фиксирует
4 рабочих потока
внутри теста, поэтому масштабирование этим тестом не измеряется. Для ALL в таблице
используется конфигурация `total_workers = ranks × threads`.

### Команды сборки

```bash
git submodule update --init --recursive
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build --config Release --target ppc_func_tests ppc_perf_tests --parallel
```

Из-за того, что реализации находятся в разных ветках, фактические запуски выполнялись
в соответствующих build-каталогах: `build_omp_verify`, `build_tbb_verify_only`,
`build_stl_verify`, `build_all_verify` и временном SEQ worktree `build_report`.

---

## 4. Сводка корректности

| Реализация |                Команда/фильтр                |    Результат    |
| ---------- | -------------------------------------------- | --------------- |
|      SEQ   | `GusevDoubleSortEvenOddBatcherSEQ.*`         | 21 tests passed |
|      OMP   | `*gusev_d_double_sort_even_odd_batcher_omp*` | 21 tests passed |
|      TBB   | `*gusev_d_double_sort_even_odd_batcher_tbb*` | 22 tests passed |
|      STL   | `*gusev_d_double_sort_even_odd_batcher_stl*` | 23 tests passed |
|      ALL   | `*gusev_d_double_sort_even_odd_batcher_all*` | 23 tests passed |

Все реализации сравнивают результат с `std::ranges::sort`. Дополнительно проверяются
крайние случаи: пустой массив, один элемент, отсортированный вход, обратный порядок,
повторы, случайные размеры, экстремальные значения, вход меньше числа потоков и
сохранение входного снимка после `PreProcessing`.

---

## 5. Агрегированные результаты

Эта таблица сводит фактические запуски от 04.06.2026. Размеры указаны явно, потому что
наборы тестов производительности в отдельных ветках исторически отличаются.

| Реализация | Режим    | Размер  | workers | time (s) | speedup | efficiency | Примечание                |
| ---------- | ---      | ---:    | ---:    | ---:     | ---:    | ---:       | --- |
|    SEQ     | task_run | 3 000   | 1       | 0.002658 | 1.00    | 1.00       | случайный вход            |
|    OMP     | Run      | 8 192   | 1       | 0.049106 | 1.00    | 1.00       | случайный вход |
|    OMP     | Run      | 8 192   | 2       | 0.041313 | 1.19    | 0.59       | случайный вход |
|    OMP     | Run      | 8 192   | 4       | 0.056564 | 0.87    | 0.22       | случайный вход |
|    OMP     | Run      | 8 192   | 8       | 0.051755 | 0.95    | 0.12       | случайный вход |
|    TBB     | Run      | 8 192   | 1       | 0.000536 | 1.00    | 1.00       | случайный вход |
|    TBB     | Run      | 8 192   | 2       | 0.085974 | 0.0062  | 0.0031     | случайный вход |
|    TBB     | Run      | 8 192   | 4       | 0.129266 | 0.0041  | 0.0010     | случайный вход |
|    TBB     | Run      | 8 192   | 8       | 0.177895 | 0.0030  | 0.0004     | случайный вход            |
|    STL     | task_run | 32 768  | 4       | 0.002519 | 1.00    | 0.25       | зафиксировано perf-тестом |
|    ALL     | task_run | 32 768  | 4       | 0.004872 | 1.00    | 0.25       | 1 rank × 4 threads        |
|    ALL     | task_run | 32 768  | 8       | 0.004526 | 1.08    | 0.13       | 2 ranks × 4 threads       |

Отдельные значения `pipeline` и дополнительные наборы входов приведены в локальных
отчетах.

---

## 6. Интерпретация различий

SEQ является контрольной реализацией и показывает время без расходов на потоки,
планировщик или MPI. Она важна для проверки корректности, но не масштабируется.

OMP на данном размере входа не показывает линейного ускорения. Два рабочих потока дали
небольшой выигрыш, а 4 и 8 потоков уже проигрывают из-за накладных расходов OpenMP и
сравнительно малого размера массива.

TBB в текущем наборе тестов производительности оказался чувствителен к гранулярности задач. При одном worker-е
ветви с порогами часто выполняются последовательно и быстро; при увеличении числа
рабочих потоков среда выполнения TBB начинает платить за планирование задач, что видно по таблице.

STL проходит стандартные `pipeline` и `task_run`, но его perf-тест фиксирует 4 рабочих
потока через guard. Поэтому локальный отчет честно не делает вывод о масштабируемости STL по
этому стандартному тесту. Главная инженерная часть STL-версии - корректный запуск,
`join()` и передача исключений из рабочих потоков.

ALL добавляет MPI-уровень. На `mpiexec -n 2` случайный `task_run` стал немного быстрее,
но эффективность низкая: часть времени уходит на `Scatterv`, `Gatherv`, `Bcast` и
финальное слияние на rank 0.

---

## 7. Репродуцируемость

Запуск функциональных тестов:

```bash
PPC_NUM_THREADS=4 ./build_stl_verify/bin/ppc_func_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_stl*"
PPC_NUM_THREADS=4 ./build_omp_verify/bin/ppc_func_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_omp*"
PPC_NUM_THREADS=4 ./build_tbb_verify_only/bin/ppc_func_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_tbb*"
PPC_NUM_THREADS=4 PPC_NUM_PROC=1 ./build_all_verify/bin/ppc_func_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_all*"
```

Запуск тестов производительности:

```bash
PPC_NUM_THREADS=4 ./build_stl_verify/bin/ppc_perf_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_stl*"
PPC_NUM_THREADS=4 ./build_omp_verify/bin/ppc_perf_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_omp*"
PPC_NUM_THREADS=4 ./build_tbb_verify_only/bin/ppc_perf_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_tbb*"
PPC_NUM_THREADS=4 PPC_NUM_PROC=1 ./build_all_verify/bin/ppc_perf_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_all*"
PPC_NUM_THREADS=4 mpiexec -n 2 ./build_all_verify/bin/ppc_perf_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher_all*"
```

SEQ запускался в отдельном worktree ветки `gusev_d_double_sort_even_odd_batcher`:

```bash
./build_report/bin/ppc_func_tests.exe --gtest_filter="GusevDoubleSortEvenOddBatcherSEQ.*"
./build_report/bin/ppc_perf_tests.exe --gtest_filter="*gusev_d_double_sort_even_odd_batcher*"
```

Дополнительной стабилизации CPU вроде отключения частотного скейлинга не выполнялось,
поэтому малые времена следует рассматривать как локальные измерения, а не как
лабораторно изолированный бенчмарк.

---

## 8. Заключение

Все пять реализаций корректно сортируют `double` и проходят свои функциональные и
тесты производительности. Самая надежная часть работы - общий алгоритмический инвариант:
преобразование `double` в сортируемый ключ перед radix sort.

На измеренных размерах входа параллельные версии не всегда ускоряются, потому что задача
быстро становится ограниченной накладными расходами на потоки, планировщик TBB, MPI-обмен
и слияние блоков. Это важный вывод: корректная параллельная реализация сама по себе не
гарантирует ускорение на малом тестовом размере.

---

## 9. Источники

1. `report-instructions.md` в корне репозитория.
2. Документация курса и структура `tasks/example_threads`.
3. OpenMP API Specification.
4. oneTBB specification and API reference.
5. MPI Forum: MPI Standard.
6. Microsoft Learn: C++ Standard Library threads.

---

## 10. Приложение

Локальные отчеты:

- `seq/report.md`
- `omp/report.md`
- `tbb/report.md`
- `stl/report.md`
- `all/report.md`

Схема файлов в ветке отчетов:

```text
tasks/gusev_d_double_sort_even_odd_batcher_reports/
  report.md
  info.json
  settings.json
  seq/report.md
  omp/report.md
  tbb/report.md
  stl/report.md
  all/report.md
```
