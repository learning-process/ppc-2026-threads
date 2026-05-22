# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников)

- Student: Назарова Ксения Олеговна, group 3823Б1ПР3
- Variant: 9
- Task directory: `tasks/Nazarova_K_rad_sort_batcher_metod`
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

В работе реализовано численное интегрирование функции нескольких переменных методом средних прямоугольников. Область интегрирования задаётся нижними и верхними границами по каждой координате, а также числом разбиений по каждой оси. Метод хорошо подходит для сравнения разных моделей параллелизма, потому что вычисление значения функции в каждой ячейке сетки независимо от остальных ячеек, а итоговый результат получается суммированием частичных вкладов.

## 2. Единая постановка задачи

Входные данные:

- `function` — интегрируемая функция вида `double f(const std::vector<double>&)`;
- `lower_bounds` — нижние границы интегрирования;
- `upper_bounds` — верхние границы интегрирования;
- `steps` — число прямоугольников по каждой координате.

Выходные данные: приближённое значение многомерного интеграла.

Для каждой размерности `i` вычисляется шаг

```text
h_i = (upper_bounds[i] - lower_bounds[i]) / steps[i]
```

Затем перебираются все ячейки декартовой сетки. Вклад ячейки считается в её середине, после чего сумма значений функции умножается на объём одной ячейки:

```text
I ≈ (sum f(x_center)) * product(h_i)
```

Корректный вход должен содержать непустую размерность, заданную функцию, согласованные размеры векторов, конечные значения границ, положительное число шагов и условие `lower_bounds[i] <= upper_bounds[i]` для каждой координаты. В параллельных версиях дополнительно проверяется отсутствие переполнения при вычислении общего числа ячеек.

## 3. Единая методика эксперимента

Эксперименты выполнялись в уже собранном каталоге `build` на следующем окружении:

| Параметр | Значение |
|---|---|
| OS | Linux 6.6.114.1-microsoft-standard-WSL2 x86_64 |
| CPU | 12th Gen Intel(R) Core(TM) i5-1235U |
| Логические CPU | 12 |
| Ядра / потоки | 6 cores, 2 threads per core |
| RAM | 7.6 GiB |
| Compiler | `c++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0` |
| CMake | 3.28.3 |
| Основная потоковая конфигурация | `PPC_NUM_THREADS=4`, `OMP_NUM_THREADS=4` |
| MPI-конфигурация для `all` | `mpirun -np 2`, дополнительно 4 worker-потока на процесс |

Функциональные тесты используют общий набор из 7 случаев: константная функция в 2D, линейная функция в 1D, произведение координат в 2D, сумма координат в 3D, квадрат в 1D, тригонометрическая функция в 2D и сдвинутое произведение в 3D. Performance-тест использует интеграл на области `[0,1] × [0,2] × [0,3]` с сеткой `120 × 120 × 120`; интегрируемая функция `x + 2y + 3z`, ожидаемый результат равен `42.0`.

Ускорение считается относительно соответствующего времени `seq`:

```text
speedup = T_seq / T_backend
efficiency = speedup / workers
```

Для `omp`, `tbb` и `stl` число работников в таблице принято равным 4. Для `all` в однопроцессном performance-запуске указано `1 × 4`, а для отдельного MPI-запуска — `2 × 4`.

## 4. Сводка корректности

Функциональный запуск для `seq`, `omp`, `tbb` и `stl`:

```bash
PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 ./build/bin/ppc_func_tests \
  --gtest_filter='RectangleIntegrationTests/NazarovaKCalcIntegRectanglesRunFuncTests.*'
```

Результат: 35 тестов из набора были запущены; 28 тестов для `seq/omp/tbb/stl` прошли, 7 тестов для `all` были пропущены, потому что `all` и MPI-задачи должны запускаться под `mpirun`.

Функциональный запуск для `all`:

```bash
mpirun --allow-run-as-root -np 2 ./build/bin/ppc_func_tests \
  --gtest_filter='RectangleIntegrationTests/NazarovaKCalcIntegRectanglesRunFuncTests.*all_enabled*'
```

Результат: все 7 тестов `all` прошли.

Таким образом, все реализации проверены на одном наборе аналитически заданных интегралов и сравниваются с теми же ожидаемыми значениями.

## 5. Агрегированные результаты

Основной performance-запуск:

```bash
PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/NazarovaKCalcIntegRectanglesRunPerfTests.*'
```

| Backend | Mode | Workers | Time, s | Speedup vs SEQ | Efficiency | Notes |
|---|---|---:|---:|---:|---:|---|
| SEQ | pipeline | 1 | 0.7317971416 | 1.00 | 1.00 | baseline |
| OMP | pipeline | 4 | 1.7958475872 | 0.41 | 0.10 | `parallel for`, `reduction` |
| TBB | pipeline | 4 | 1.1344363388 | 0.65 | 0.16 | `parallel_reduce` |
| STL | pipeline | 4 | 2.0161764572 | 0.36 | 0.09 | `std::async` chunks |
| ALL | pipeline | 1 × 4 | 1.9185 | 0.38 | 0.10 | MPI size 1 + TBB |
| SEQ | task_run | 1 | 0.6949633072 | 1.00 | 1.00 | baseline |
| OMP | task_run | 4 | 1.7755917146 | 0.39 | 0.10 | `parallel for`, `reduction` |
| TBB | task_run | 4 | 1.0996614686 | 0.63 | 0.16 | `parallel_reduce` |
| STL | task_run | 4 | 1.9265005476 | 0.36 | 0.09 | `std::async` chunks |
| ALL | task_run | 1 × 4 | 1.9017 | 0.37 | 0.09 | MPI size 1 + TBB |

Отдельный запуск `all` под MPI:

```bash
PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 mpirun --allow-run-as-root -np 2 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/NazarovaKCalcIntegRectanglesRunPerfTests.*all_enabled*'
```

| Backend | Mode | Ranks | Threads per rank | Total workers | Time, s | Speedup vs SEQ | Efficiency |
|---|---|---:|---:|---:|---:|---:|---:|
| ALL | pipeline | 2 | 4 | 8 | 1.0314585568 | 0.71 | 0.09 |
| ALL | task_run | 2 | 4 | 8 | 1.0240943998 | 0.68 | 0.08 |

## 6. Интерпретация различий

Последовательная версия оказалась быстрее чистых потоковых запусков на выбранном размере. Это не означает ошибку в параллельных реализациях: вычисление значения функции в тесте очень лёгкое, поэтому накладные расходы на создание задач, разбиение диапазона, редукцию и работу runtime оказываются сопоставимы с полезной работой.

OpenMP-версия использует статическое разбиение линейного диапазона ячеек и редукцию суммы. Такая схема проста и детерминирована по структуре работы, но при лёгком теле цикла цена параллельной области и редукции снижает итоговую эффективность.

TBB-версия использует `parallel_reduce`. Она удобна для суммирования независимых вкладов, но на данном тесте runtime oneTBB и создание локальных диапазонов дают больший overhead, чем выигрыш от параллельного обхода.

STL-версия вручную делит диапазон между асинхронными задачами и собирает локальные суммы через `future::get`. Среди чисто потоковых backend-ов она показала лучший результат в этом запуске, но всё равно не превзошла последовательный baseline из-за цены запуска асинхронных задач.

Гибридная версия `all` распределяет линейный диапазон между MPI-процессами, внутри каждого процесса применяет TBB и затем собирает сумму через `MPI_Allreduce`. В однопроцессном режиме она ведёт себя близко к TBB с дополнительной MPI-обвязкой. При `2 × 4` время заметно уменьшается относительно `1 × 4`, но эффективность по общему числу работников остаётся низкой из-за коммуникационных и runtime-издержек.

## 7. Репродуцируемость

Сборка:

```bash
cmake -S . -B build \
  -D USE_FUNC_TESTS=ON \
  -D USE_PERF_TESTS=ON \
  -D CMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Функциональные тесты:

```bash
PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 ./build/bin/ppc_func_tests \
  --gtest_filter='RectangleIntegrationTests/NazarovaKCalcIntegRectanglesRunFuncTests.*'

mpirun --allow-run-as-root -np 2 ./build/bin/ppc_func_tests \
  --gtest_filter='RectangleIntegrationTests/NazarovaKCalcIntegRectanglesRunFuncTests.*all_enabled*'
```

Тесты производительности:

```bash
PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/NazarovaKCalcIntegRectanglesRunPerfTests.*'

PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 mpirun --allow-run-as-root -np 2 ./build/bin/ppc_perf_tests \
  --gtest_filter='RunModeTests/NazarovaKCalcIntegRectanglesRunPerfTests.*all_enabled*'
```

## 8. Заключение

Все реализации решают одну и ту же задачу и проходят функциональные тесты. На выбранном performance-тесте последовательная версия остаётся самой быстрой, потому что интегрируемая функция слишком дешева относительно накладных расходов параллельных runtime. Гибридная версия показывает улучшение при запуске на двух MPI-процессах по сравнению с однопроцессным `all`, однако абсолютное ускорение относительно `seq` всё ещё меньше единицы.

Для получения более выраженного ускорения следует использовать более дорогую интегрируемую функцию, больший размер сетки или серию повторов с более строгой стабилизацией окружения. Текущий отчёт фиксирует фактические результаты на доступной машине и не переносит их автоматически на другие конфигурации.

## 9. Источники

- Документация курса и `report-instructions.md` из репозитория.
- OpenMP Application Programming Interface Specification.
- oneTBB documentation: `parallel_reduce`, `blocked_range`.
- MPI standard documentation: `MPI_Comm_rank`, `MPI_Comm_size`, `MPI_Allreduce`.
- C++ reference documentation for `std::async`, `std::future` and threading primitives.
