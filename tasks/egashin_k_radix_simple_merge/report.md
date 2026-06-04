# Поразрядная сортировка для вещественных чисел double с простым слиянием

- Student: Егашин Кирилл Олегович, group 3823Б1ФИ2
- Variant: 19
- Technology reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

Работа сравнивает пять реализаций одной задачи: SEQ, OMP, TBB, STL и ALL. Задача хорошо подходит для сравнения
моделей параллелизма, потому что массив можно разделить на независимые части, отсортировать их отдельно и
затем объединить простым слиянием.

## 2. Постановка задачи

Вход: массив `std::vector<double>`.

Выход: отсортированный по возрастанию массив той же длины.

Критерий корректности: выход должен совпадать с результатом `std::ranges::sort` для тех же входных данных.

## 3. Единая методика

Для всех реализаций используется один базовый подход:

1. Преобразовать `double` в сортируемый 64-битный ключ.
2. Выполнить LSD radix sort по байтам ключа.
3. Для параллельных версий разделить вход на блоки.
4. Отсортировать блоки независимо.
5. Объединить отсортированные блоки простым merge.

SEQ сортирует весь массив одним блоком. OMP, TBB и STL сортируют блоки в потоках. ALL дополнительно делит
данные между MPI рангами и внутри каждого ранга использует OpenMP.

## 4. Сводка реализации

Основные файлы:

- `common/include/common.hpp`
- `common/include/radix_utils.hpp`
- `seq/src/ops_seq.cpp`
- `omp/src/ops_omp.cpp`
- `tbb/src/ops_tbb.cpp`
- `stl/src/ops_stl.cpp`
- `all/src/ops_all.cpp`
- `tests/functional/main.cpp`
- `tests/performance/main.cpp`

Общий helper `radix_utils.hpp` используется новыми параллельными версиями для сортировки диапазона и простого
слияния. SEQ оставлена без рефакторинга, чтобы не менять уже готовую базовую реализацию.

## 5. Корректность

Функциональные тесты покрывают:

- пустой массив
- один элемент
- уже отсортированный массив
- обратный порядок
- смешанные знаки
- повторы
- разные масштабы значений
- положительную и отрицательную бесконечность

Performance test проверяет, что результат отсортирован, через `std::ranges::is_sorted`.

## 6. Экспериментальная среда

Окружение:

- OS: Ubuntu 24.04
- Compiler: `g++-14`
- Build type: `Release`
- Build system: CMake + Ninja
- Threads: `PPC_NUM_THREADS=2`
- Processes: `PPC_NUM_PROC=2`
- Metric: `task_run`

## 7. Команды воспроизведения

```bash
cmake -S . -B build -G Ninja -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel
scripts/run_tests.py --running-type=threads --counts 1 2 4 --build-dir build
PPC_NUM_THREADS=2 scripts/run_tests.py --running-type=processes --counts 1 2 4 --build-dir build
PPC_NUM_THREADS=2 PPC_NUM_PROC=2 scripts/run_tests.py --running-type=performance --build-dir build
```

## 8. Агрегированные результаты

| Backend | Configuration     | Time, s  | Speedup vs seq | Efficiency |
| ------- | ----------------- | -------- | -------------- | ---------- |
| seq     | 1 worker          | 0.035714 | 1.00           | N/A        |
| omp     | 2 threads         | 0.024104 | 1.48           | 0.74       |
| tbb     | 2 workers         | 0.023022 | 1.55           | 0.78       |
| stl     | 2 threads         | 0.024620 | 1.45           | 0.73       |
| all     | 2 ранга, 2 потока | 0.030749 | 1.16           | 0.29       |

## 9. Интерпретация

SEQ задает baseline и проверочный результат.

OMP ожидаемо должен давать выигрыш на больших массивах, потому что сортировка блоков независима.

TBB использует ту же блочную схему, но планирование выполняет runtime TBB.

STL дает ручной контроль над потоками, но платит за создание и синхронизацию `std::thread`.

ALL добавляет процессный уровень. На больших входах это может помочь, но на малых входах коммуникации и
финальное слияние на ранге 0 могут перекрыть выигрыш.

## 10. Источники

1. Документация курса: `docs/common_information/report.rst`
2. Документация курса: `docs/common_information/threading_tasks.rst`
3. OpenMP specification: <https://www.openmp.org/specifications/>
4. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
5. MPI standard: <https://www.mpi-forum.org/docs/>
