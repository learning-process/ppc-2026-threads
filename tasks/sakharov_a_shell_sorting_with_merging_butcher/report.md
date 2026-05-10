# Сортировка Шелла с четно-нечетным слиянием Бэтчера

- **Student:** Сахаров Александр Владимирович, группа 3823Б1ФИ3
- **Variant:** 16
- **Local reports:** `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

Цель работы — реализовать и сравнить несколько версий сортировки массива целых чисел: последовательную, OpenMP, oneTBB,
STL threads и гибридную ALL-версию. Задача хорошо подходит для сравнения моделей параллелизма, потому что сортировка
независимых частей массива и слияния на одном уровне дерева могут выполняться параллельно.

## 2. Единая постановка задачи

На вход всем версиям подается `std::vector<int>`. На выходе должен быть `std::vector<int>` с теми же элементами в
порядке неубывания. Проверяются пустой массив, один элемент, отрицательные значения, повторяющиеся элементы, четные и
нечетные размеры входа.

Контракт реализации:

- входной тип: `std::vector<int>`;
- выходной тип: `std::vector<int>`;
- базовый класс задачи: `ppc::task::BaseTask`;
- эталон для проверки: `std::ranges::sort`;
- критерий корректности: полное совпадение результата с эталоном.

## 3. Единая методика эксперимента

Все реализации используют один набор функциональных тестов и один набор perf-тестов. В perf-тестах используется массив
из `400000` целых чисел. Вход генерируется детерминированно, чтобы все backend-ы работали с одинаковыми данными.

Экспериментальная среда:

- **CPU:** Intel Core i5-12400F, 6 cores / 12 threads, 2.50 GHz.
- **RAM:** 32 ГБ.
- **OS:** Windows 10 + WSL, Ubuntu 24.04.3 LTS.
- **IDE:** Visual Studio Code.
- **Compiler:** `g++ 13.3.0`.
- **Build system:** `CMake 3.28.3`.
- **VCS:** Git.
- **Build type:** Release-сборка проекта.
- **Perf modes:** `pipeline` и `task_run`.

Для OpenMP, oneTBB и STL workers равны значению `PPC_NUM_THREADS`. Для ALL workers задаются как `ranks x threads`; при
расчете эффективности используется произведение этих величин. В таблице приведены средние значения для конфигураций `2`,
`4`, `8` потоков и `2 x 2`, `2 x 4`, `2 x 8` для ALL. Для конфигураций с `2` и `8` потоками среднее рассчитано по трем
запускам, для конфигурации `4` потоков — по двум запускам.

Формулы:

```text
speedup = time_seq / time_parallel
efficiency = speedup / workers
```

## 4. Сводка корректности

Функциональные тесты находятся в `tests/functional/main.cpp`. Для каждого входа ожидаемый результат строится стандартной
сортировкой, затем сравнивается с результатом выбранной реализации.

Проверяются:

- массив с четным числом элементов;
- массив с отрицательными числами и нечетным размером;
- массив с повторяющимися значениями;
- массив из одного элемента;
- пустой массив.

По результатам локального запуска все версии прошли функциональные тесты: SEQ, OpenMP, oneTBB, STL и ALL получили 5/5
успешных тестов.

## 5. Агрегированные результаты

В таблице приведены средние результаты по трем запускам perf-тестов.

| Реализация | Workers | Pipeline, с | Task run, с | Speedup pipeline | Speedup task_run | Efficiency task_run |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SEQ | 1 | 0.0448361238 | 0.0451068401 | 1.00 | 1.00 | 1.00 |
| OpenMP | 2 | 0.0141551325 | 0.0087188449 | 3.17 | 5.17 | 2.59 |
| OpenMP | 4 | 0.0106660785 | 0.0061448581 | 4.20 | 7.34 | 1.84 |
| OpenMP | 8 | 0.0172132381 | 0.0062423649 | 2.60 | 7.23 | 0.90 |
| oneTBB | 2 | 0.0118986448 | 0.0120335897 | 3.77 | 3.75 | 1.87 |
| oneTBB | 4 | 0.0078377247 | 0.0078923226 | 5.72 | 5.72 | 1.43 |
| oneTBB | 8 | 0.0066771825 | 0.0064432303 | 6.71 | 7.00 | 0.88 |
| STL | 2 | 0.0256135146 | 0.0233078480 | 1.75 | 1.94 | 0.97 |
| STL | 4 | 0.0140028239 | 0.0143011809 | 3.20 | 3.15 | 0.79 |
| STL | 8 | 0.0106910864 | 0.0100878557 | 4.19 | 4.47 | 0.56 |
| ALL | 2 x 2 | 0.0100328291 | 0.0088704593 | 4.47 | 5.09 | 1.27 |
| ALL | 2 x 4 | 0.0094135525 | 0.0090741707 | 4.76 | 4.97 | 0.62 |
| ALL | 2 x 8 | 0.0101951866 | 0.0095393695 | 4.40 | 4.73 | 0.30 |

## 6. Интерпретация различий

SEQ является baseline-версией и не использует параллелизм. OpenMP показывает лучший результат в режиме `task_run`:
основной вычислительный этап хорошо раскладывается на независимую сортировку чанков и параллельные проходы слияния.
oneTBB показывает лучшее среднее время в режиме `pipeline`, что связано с эффективным планированием независимых задач
runtime-библиотекой.

STL-версия ускоряет вычисления относительно SEQ, но проигрывает OpenMP и oneTBB из-за ручного создания и синхронизации
потоков. ALL-версия также быстрее SEQ, но ее эффективность падает при росте общего числа workers из-за MPI-обменов,
сборки данных на rank 0 и финального слияния собранных частей. Значения эффективности выше `1` на некоторых малых
конфигурациях объясняются тем, что параллельные версии используют другую организацию сортировки чанков по сравнению с
SEQ baseline, а также естественной вариативностью коротких замеров.

Отдельно это заметно на `2` и `4` потоках: OpenMP и oneTBB получают выигрыш не только от параллельного выполнения, но и
от перехода к схеме `std::sort` по чанкам с деревом `std::merge`. Последовательная версия остается корректным baseline
для задачи, но она не является однопоточным вариантом той же самой chunk-based схемы. Поэтому эффективность в таблице
является относительной метрикой для сравнения реализаций в рамках работы, а не чистой аппаратной эффективностью потоков.

## 7. Репродуцируемость

Сборка:

```bash
cmake --build build --target ppc_func_tests ppc_perf_tests -j
```

Функциональные тесты:

```bash
./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_seq*'

PPC_NUM_THREADS=8 ./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_omp*'

PPC_NUM_THREADS=8 ./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_tbb*'

PPC_NUM_THREADS=8 ./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_stl*'

PPC_NUM_THREADS=8 mpirun -np 2 ./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_all*'
```

Perf-тесты:

```bash
./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_seq*'

PPC_NUM_THREADS=8 ./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_omp*'

PPC_NUM_THREADS=8 ./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_tbb*'

PPC_NUM_THREADS=8 ./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_stl*'

PPC_NUM_THREADS=8 mpirun -np 2 ./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_all*'
```

## 8. Заключение

Все требуемые версии задачи реализованы и прошли функциональные тесты. По усредненным измерениям oneTBB лучше всего
показала себя в режиме `pipeline`, а OpenMP — в режиме `task_run`. ALL-версия корректно использует MPI, OpenMP, oneTBB и
STL threads, но на текущем размере входа часть выигрыша теряется на коммуникации и финальном слиянии.

## 9. Источники

- Лекции Сысоева А. В. по курсу «Параллельное программирование для систем с общей памятью».
- [OpenMP API Specification](https://www.openmp.org/specifications/).
- [oneAPI Threading Building Blocks Documentation](https://uxlfoundation.github.io/oneTBB/).
- [MPI Standard Documentation](https://www.mpi-forum.org/docs/).
- cppreference: [`std::thread`](https://en.cppreference.com/w/cpp/thread/thread),
  [`std::sort`](https://en.cppreference.com/w/cpp/algorithm/sort),
  [`std::merge`](https://en.cppreference.com/w/cpp/algorithm/merge).

## 10. Приложение

Ключевой фрагмент ALL-версии:

```cpp
auto local_data = ScatterInput(input, root_data, rank);
local_data = SortLocalPart(std::move(local_data));
auto gathered_data = GatherSortedData(local_data, root_data, rank);
auto result = MergeRootData(std::move(gathered_data), root_data, rank);
return BroadcastResult(std::move(result), rank);
```

Этот фрагмент отражает основную схему гибридной реализации: MPI распределяет и собирает данные, локальная часть
сортируется внутри процесса, а готовый результат рассылается всем rank-ам.
