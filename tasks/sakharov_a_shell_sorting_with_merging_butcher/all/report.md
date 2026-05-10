# Сортировка Шелла с четно-нечетным слиянием Бэтчера — ALL

- **Student:** Сахаров Александр Владимирович, группа 3823Б1ФИ3
- **Technology:** ALL
- **Variant:** 16

## 1. Контекст

ALL-версия объединяет несколько технологий параллельного программирования: MPI, OpenMP, Intel oneTBB и STL threads.
Такая реализация нужна, чтобы показать одновременно межпроцессное распределение данных и внутрипроцессную параллельную
обработку локальных частей массива.

## 2. Постановка задачи

На вход поступает массив целых чисел `std::vector<int>`. Нужно получить массив того же типа, содержащий те же элементы в
порядке неубывания. ALL-версия должна совпадать по результату с SEQ-версией и корректно работать при запуске через
`mpirun`.

## 3. Базовый алгоритм

Базовая логика остается такой же, как у остальных параллельных версий: вход делится на части, части сортируются
независимо, затем объединяются деревом слияний. В ALL-версии это разбиение выполняется в два уровня: сначала между
MPI-процессами, затем внутри каждого процесса между потоками.

## 4. Межпроцессная схема

На уровне MPI работа строится так:

1. Нулевой rank вычисляет размеры частей массива для всех процессов.
2. `MPI_Scatter` передает каждому процессу размер его локального блока.
3. `MPI_Scatterv` распределяет сами элементы массива.
4. Каждый процесс сортирует свою локальную часть.
5. `MPI_Gatherv` собирает локально отсортированные части на rank 0.
6. Rank 0 выполняет финальное слияние собранных частей.
7. `MPI_Bcast` рассылает итоговый отсортированный массив всем процессам.

## 5. Внутрипроцессная схема

Внутри каждого MPI-процесса используются OpenMP и oneTBB:

- OpenMP сортирует локальные чанки через `#pragma omp parallel for`;
- oneTBB выполняет локальное дерево слияний через `tbb::parallel_for`;
- STL threads используются на rank 0 для финального слияния частей, собранных от MPI-процессов.

Конфигурация запуска фиксируется как `ranks x threads`. В текущих замерах использовались запуски `2 x 2`, `2 x 4` и `2 x
8`: два MPI-процесса и разное число потоков на процесс.

## 6. Детали реализации

Файлы реализации:

- `all/include/ops_all.hpp`;
- `all/src/ops_all.cpp`.

Класс `SakharovAShellButcherALL` возвращает тип задачи `ppc::task::TypeOfTask::kALL`.

Минимальный листинг из `all/src/ops_all.cpp`:

```cpp
std::vector<int> RunMpiSort(const std::vector<int> &input, int rank, int process_count) {
  const auto root_data = BuildRootData(input, rank, process_count);
  auto local_data = ScatterInput(input, root_data, rank);
  local_data = SortLocalPart(std::move(local_data));
  auto gathered_data = GatherSortedData(local_data, root_data, rank);
  auto result = MergeRootData(std::move(gathered_data), root_data, rank);
  return BroadcastResult(std::move(result), rank);
}
```

Этот фрагмент показывает полный жизненный цикл ALL-версии: распределение данных, локальную сортировку, сбор, финальное
слияние и рассылку результата.

Дополнительный фрагмент с MPI-обменом:

```cpp
MPI_Scatterv(RootInputBuffer(input, rank), RootBuffer(root_data.counts, rank),
             RootBuffer(root_data.displacements, rank), MPI_INT,
             BufferOrNull(local_data), local_size, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Gatherv(BufferOrNull(local_data), local_size, MPI_INT, BufferOrNull(gathered_data),
            RootBuffer(root_data.counts, rank), RootBuffer(root_data.displacements, rank),
            MPI_INT, 0, MPI_COMM_WORLD);
```

`MPI_Scatterv` нужен, потому что части массива могут иметь разный размер. `MPI_Gatherv` по той же причине используется
при сборе отсортированных локальных частей.

## 7. Проверка корректности

Корректность ALL-версии проверяется общими функциональными тестами под `mpirun`. Эталонный результат строится через
`std::ranges::sort`. После `MPI_Bcast` итоговый массив доступен всем MPI-процессам, поэтому проверка результата
выполняется одинаково на каждом rank-е. По результатам локального запуска все 5 функциональных тестов ALL прошли
успешно.

## 8. Экспериментальная среда

- **CPU:** Intel Core i5-12400F, 6 cores / 12 threads, 2.50 GHz.
- **RAM:** 32 ГБ.
- **OS:** Windows 10 + WSL, Ubuntu 24.04.3 LTS.
- **IDE:** Visual Studio Code.
- **Compiler:** `g++ 13.3.0`.
- **Build system:** `CMake 3.28.3`.
- **VCS:** Git.
- **Build type:** Release-сборка проекта.
- **Perf modes:** `pipeline` и `task_run`.
- **Hybrid configuration:** в таблице приведены результаты для `2 x 2`, `2 x 4` и `2 x 8`.

## 9. Результаты

В таблице приведены средние результаты perf-тестов. Для `2 x 2` и `2 x 8` среднее рассчитано по трем запускам, для `2 x
4` — по двум запускам. Эффективность считается относительно общего числа workers: `ranks * threads`.

| Режим | Ranks x Threads | Время, с | Ускорение относительно SEQ | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| pipeline | 2 x 2 | 0.0100328291 | 4.47 | 1.12 |
| task_run | 2 x 2 | 0.0088704593 | 5.09 | 1.27 |
| pipeline | 2 x 4 | 0.0094135525 | 4.76 | 0.60 |
| task_run | 2 x 4 | 0.0090741707 | 4.97 | 0.62 |
| pipeline | 2 x 8 | 0.0101951866 | 4.40 | 0.27 |
| task_run | 2 x 8 | 0.0095393695 | 4.73 | 0.30 |

Для конфигурации `2 x 2` эффективность выше `1`, потому что метрика считается относительно SEQ baseline, а ALL-версия
использует другую организацию работы: MPI-разбиение входа, локальную сортировку чанков через `std::sort`,
OpenMP-параллелизм, oneTBB-слияние и финальное STL-слияние. Поэтому ускорение отражает не только увеличение числа
workers, но и отличие вычислительной схемы от последовательной реализации. При росте числа workers до `2 x 4` и `2 x 8`
эффективность снижается из-за стоимости MPI-обмена, сборки результата и финального слияния на rank 0.

## 10. Репродуцируемость

Функциональные тесты:

```bash
PPC_NUM_THREADS=8 mpirun -np 2 ./build/bin/ppc_func_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_all*'
```

Perf-тесты:

```bash
PPC_NUM_THREADS=8 mpirun -np 2 ./build/bin/ppc_perf_tests \
  --gtest_filter='*sakharov_a_shell_sorting_with_merging_butcher_all*'
```

## 11. Выводы

ALL-версия корректно объединяет MPI, OpenMP, oneTBB и STL threads. Ускорение относительно SEQ есть, но эффективность
ниже, чем у чистых OpenMP и oneTBB, из-за стоимости распределения данных, сборки результата и финального слияния на rank
0. Такой подход лучше раскрывается на более крупных входах, где коммуникационные расходы становятся менее заметными
относительно вычислений.

## 12. Источники

- Лекции Сысоева А. В. по курсу «Параллельное программирование для систем с общей памятью».
- [MPI Standard Documentation](https://www.mpi-forum.org/docs/).
- [OpenMP API Specification](https://www.openmp.org/specifications/).
- [oneAPI Threading Building Blocks Documentation](https://uxlfoundation.github.io/oneTBB/).
- cppreference: [`std::thread`](https://en.cppreference.com/w/cpp/thread/thread),
  [`std::sort`](https://en.cppreference.com/w/cpp/algorithm/sort),
  [`std::merge`](https://en.cppreference.com/w/cpp/algorithm/merge).
