# Вычисление многомерных интегралов методом прямоугольников — ALL

- Студент: Дергунов Сергей Антонович, группа 3823Б1ПР4
- Технология: ALL
- Вариант: 9

## 1. Введение

ALL-версия реализована как гибридная MPI+OpenMP версия алгоритма.
Для распределения работы между процессами используется MPI,
а внутри каждого процесса — OpenMP для распараллеливания по потокам.

## 2. Постановка задачи

- **Вход**:
  - подынтегральная функция `std::function<double(const std::vector<double>&)>`
  - границы интегрирования `std::vector<std::pair<double, double>>`
  - количество шагов разбиения `int`
- **Условие корректности входа**:
  - функция задана корректно
  - количество шагов > 0
  - границы не пусты и для каждого измерения `left < right`
- **Ограничения**: размерность произвольная, но для тестов используется 1-3
- **Выход**: приближённое значение интеграла `double`

## 3. Базовый алгоритм (последовательный)

Алгоритм совпадает с SEQ (метод средних прямоугольников).

## 4. Схема распараллеливания

ALL-версия использует гибридную схему MPI + OpenMP:

- **MPI уровень**: общее количество точек сетки `total_points = n^dim` равномерно
  распределяется между MPI-процессами с учётом остатка. Каждый процесс получает
  непрерывный диапазон точек для обработки.
- **OpenMP уровень**: внутри каждого MPI-процесса используется
  `#pragma omp parallel for reduction(+ : local_sum)` для распараллеливания
  обработки полученного диапазона по потокам.
- **Сбор результатов**: после завершения вычислений выполняется
  `MPI_Allreduce` для суммирования частичных результатов со всех процессов.

Конфигурация для потоковой задачи: `ranks × threads = 2 × N`, где
`ranks = 2` (количество MPI-процессов), `threads = N` (количество OpenMP-потоков).

## 5. Детали реализации

- Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`
- Класс: `DergynovSIntegralsMultistepRectangleALL`
- Конвейер: `ValidationImpl`, `PreProcessingImpl`, `RunImpl`, `PostProcessingImpl`
- Ключевые вызовы: `MPI_Comm_rank`, `MPI_Comm_size`, `MPI_Allreduce`, OpenMP directives

## 6. Экспериментальная среда

Сборка:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build --parallel
```

Запуск:

```bash
mpirun --allow-run-as-root -n 2 ./build/bin/ppc_func_tests --gtest_filter="*all*"

export PPC_NUM_THREADS=2
mpirun --allow-run-as-root -n 2 ./build/bin/ppc_perf_tests --gtest_filter="*all*"

export PPC_NUM_THREADS=4
mpirun --allow-run-as-root -n 2 ./build/bin/ppc_perf_tests --gtest_filter="*all*"
```

## 7. Результаты

### 7.1 Корректность

ALL-версия проходит все функциональные тесты (19 тестов).
Результаты совпадают с baseline SEQ с точностью до 1e-6.

### 7.2 Производительность

Определения метрик:

- `workers` — число потоков
- `time` — wall-clock время, секунды
- `speedup = T_seq / T_mode`
- `efficiency = speedup / workers * 100%`

Результаты:

| mode | workers | time (task_run), s | time (pipeline), s | speedup | efficiency, % |
|------|--------:|-------------------:|-------------------:|--------:|--------------:|
| seq  |       1 |            0.02584 |            0.02783 |    1.00 |           N/A |
| all  |       2 |            0.03955 |            0.03427 |    0.65 |         32.5% |
| all  |       8 |            0.03676 |            0.04662 |    0.70 |          8.8% |

Комментарий: ALL-версия показывает ускорение до 0.70x при 8 workers.
Накладные расходы на MPI-коммуникацию и OpenMP-синхронизацию ограничивают
масштабируемость для данной вычислительной задачи.

## 8. Выводы

ALL-версия демонстрирует гибридную модель параллелизма MPI+OpenMP.
Корректность подтверждена функциональными тестами. Ускорение ограничено
накладными расходами на коммуникацию и синхронизацию.

## 9. Источники

1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. OpenMP specification: <https://www.openmp.org/specifications/>
3. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
4. C++ reference (std::thread): <https://en.cppreference.com/w/cpp/thread/thread>
5. Course report requirements: `docs/common_information/report.rst`
