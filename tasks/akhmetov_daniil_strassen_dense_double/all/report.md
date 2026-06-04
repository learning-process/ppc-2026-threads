# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Студент: Ахметов Даниил Данисович, группа 3823Б1ПР2
- Технология: ALL
- Вариант: 3

## 1. Введение

В версии `ALL` для семестра потоков нужно показать использование
всех технологий параллелизма в одной реализации: OpenMP, `std::thread`,
oneTBB и MPI.

Основное умножение выполняется по OMP-схеме; после вычисления `C`
запускаются демонстрационные блоки для каждой технологии.
`ALL` оценивается как отдельный backend, а не как самый быстрый вариант.

Общий контекст - в [report.md](../report.md).

## 2. Постановка задачи

Та же постановка, что и в `SEQ`: `C = A * B`, допуск `1e-7`.
Дополнительно - демонстрация MPI и всех thread-технологий
в одном `RunImpl()`.

## 3. Описание алгоритма (базового/последовательного)

Вычислительная часть повторяет `OMP`: итеративный Штрассен,
`kThreshold = 64`, padding, OpenMP внутри стека.
Подробнее - [omp/report.md](../omp/report.md).

## 4. Схема распараллеливания

### MPI (межпроцессный уровень)

- `MPI_Comm_rank` - определение rank-а;
- каждый rank выполняет полное умножение (без распределения `M1..M7`);
- `MPI_Barrier(MPI_COMM_WORLD)` в конце `RunImpl()`.

### Потоки (внутрипроцессный уровень)

**Основное вычисление (OpenMP):**

- как в `OMP`: `parallel for` в умножении, split/merge, padding;
- `if (size >= kParallelThreshold)`, `kParallelThreshold = 256`.

**Демонстрация технологий (после вычисления C):**

1. OpenMP - `#pragma omp parallel` с atomic-счётчиком (`rank == 0`);
2. STL - `num_threads` потоков `std::thread` + `atomic`;
3. oneTBB - `tbb::parallel_for` по `[0, num_threads)`.

Демонстрационные блоки не изменяют матрицу `C`.

## 5. Детали реализации

**Файлы:**

- `all/include/ops_all.hpp`
- `all/src/ops_all.cpp`

**Пайплайн:**

- `ValidationImpl()` -> `PreProcessingImpl()` -> `RunImpl()` ->
  `PostProcessingImpl()`;
- `RunImpl()`: OMP-Штрассен -> демо OpenMP/STL/TBB -> `MPI_Barrier`.

**Особенности:**

- число потоков - `ppc::util::GetNumThreads()`;
- для STL используется `threads.at(i)`;
- результат умножения не зависит от демонстрационных блоков.

## 6. Экспериментальная установка

- CPU: AMD Ryzen 7 2700, 8 ядер / 16 потоков;
- RAM: 16 GB;
- OS: Windows 10 x64;
- IDE: Visual Studio Code;
- compiler: MSVC (сборка через CMake);
- build type: `Release`;
- `PPC_NUM_THREADS=4` (1 rank) или `PPC_NUM_THREADS=2`, `mpiexec -n 2`.

**Команды:**

```powershell
cd build\bin
$env:PPC_NUM_THREADS = "4"
.\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_all_*

$env:PPC_NUM_THREADS = "2"
mpiexec -n 2 .\ppc_func_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double*_all_*

$env:PPC_NUM_THREADS = "4"
.\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double_all_enabled*

$env:PPC_NUM_THREADS = "2"
mpiexec -n 2 .\ppc_perf_tests --gtest_filter=*akhmetov_daniil_strassen_dense_double_all_enabled*
```

## 7. Результаты и обсуждение

### 7.1 Корректность

- размеры `64`, `128`, `256`;
- сравнение с наивным эталоном;
- локально все тесты пройдены на `mpiexec -n 2`.

### 7.2 Производительность

| Режим | Число потоков | Время, с | Ускорение | Эффективность |
| --- | ---: | ---: | ---: | ---: |
| seq | 1 | 0.938515 | 1.00 | N/A |
| all (1 rank) | 4 | 1.662727 | 0.56 | 14.1% |
| all (2 ranks × 2 threads) | 4 | 2.187852 | 0.43 | 21.5% |

`ALL` близок к `OMP` по времени - общая вычислительная часть.
На `mpiexec -n 2` медленнее, так как каждый rank выполняет
полное умножение плюс демо-блоки и barrier.

## 8. Заключение

`ALL` корректно проходит функциональные тесты и демонстрирует
все требуемые технологии потоков и MPI. По производительности
уступает `TBB`/`STL`; его следует оценивать отдельно
с учётом конфигурации `ranks × threads`.

## 9. Источники

1. [MPI Forum Documentation](https://www.mpi-forum.org/docs/)
2. [OpenMP Specifications](https://www.openmp.org/specifications/)
3. [oneTBB Documentation](https://uxlfoundation.github.io/oneTBB/)
4. [Parallel Programming Course](https://learning-process.github.io/parallel_programming_course/ru/common_information/threading_tasks.html)

## Приложение (опционально)

```cpp
// all/src/ops_all.cpp - фрагмент после вычисления C
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if (rank == 0) { /* OpenMP demo */ }
{ /* std::thread demo */ }
{ /* tbb::parallel_for demo */ }
MPI_Barrier(MPI_COMM_WORLD);
```
