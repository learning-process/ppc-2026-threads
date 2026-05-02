# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – строковый (CRS)

- Студент: Ермаков Алексей Викторович
- Группа: 3823Б1ПР3
- Вариант: 6
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

В работе исследуется задача умножения двух разреженных комплексных матриц,
хранящихся в формате CRS. Эта задача подходит для сравнения разных моделей
параллелизма, потому что вычислительное ядро у всех реализаций одно и то же,
а различаются именно способы организации исполнения: последовательный baseline,
OpenMP, oneTBB, `std::thread` и гибридная схема `MPI + OMP`.

## 2. Единая постановка задачи

Вход:

- две матрицы `A` и `B` в формате CRS;
- каждая матрица задается полями `rows`, `cols`, `values`, `col_index`, `row_ptr`.

Выход:

- матрица `C = A * B` в формате CRS.

Ограничения и критерий корректности:

- `A.cols == B.rows`;
- `row_ptr.size() == rows + 1`;
- `values.size() == col_index.size()`;
- `row_ptr.front() == 0`;
- `row_ptr.back() == nnz`;
- все индексы `col_index[k]` лежат в диапазоне `[0, cols)`;
- корректным считается результат, совпадающий с плотным эталоном умножения.

## 3. Единая методика эксперимента

Окружение:

- CPU: AMD Ryzen 5 1600
- RAM: 8 GB
- OS: Windows 11
- compiler: MSVC
- build type: `Release`
- MPI: MS-MPI 10.0

Переменные окружения:

- `OMP_NUM_THREADS` использовалась для чистой `OMP`-версии;
- `PPC_NUM_THREADS` использовалась для `TBB`, `STL` и для OpenMP-части `ALL`;
- число MPI-процессов для `ALL` задавалось через `mpiexec -n`.

Источник данных:

- functional-тесты используют фиксированный пример `3x3` и детерминированно
  сгенерированные матрицы размеров `10x10`, `20x20`, `30x30`;
- performance-тест использует детерминированные CRS-матрицы размера
  `15000 x 15000` с плотностью `0.001`.

Размеры задач:

- для корректности: `3`, `10`, `20`, `30`;
- для производительности: `15000 x 15000`.

Метрики:

- `time` - медиана `task_run`;
- `speedup = T_seq / T_backend`;
- для `OMP`, `TBB`, `STL`: `efficiency = speedup / threads`;
- для `ALL`: `efficiency = speedup / total_workers`,
  где `total_workers = ranks * threads_per_rank`.

Серия замеров:

- для каждой конфигурации выполнялся `1` прогревочный запуск;
- после прогрева выполнялось `5` измеренных запусков;
- в таблицы вынесена медиана этих `5` значений `task_run`.

## 4. Сводка корректности

Все backend-ы сравнивались с последовательной версией `SEQ`, которая играет роль
baseline и по корректности, и по времени. Проверка выполнялась общей тестовой
инфраструктурой:

- functional-тесты из `tests/functional/main.cpp`;
- performance-тесты из `tests/performance/main.cpp`.

Использованные функциональные сценарии:

- фиксированный пример `SmallFixed` размером `3x3`;
- детерминированные разреженные матрицы `10x10`, `20x20`, `30x30`;
- сравнение с плотным эталоном `DenseMul`.

Ограничения применимости:

- `SEQ` используется как эталон;
- `OMP`, `TBB`, `STL` сравниваются по числу потоков;
- `ALL` сравнивается по конфигурации `ranks × threads`, причем в функциональном
  режиме итоговый результат дополнительно рассылается всем rank-ам.

## 5. Агрегированные результаты

Последовательный baseline:

- `T_seq = 0.1341834602 c`.

Общая таблица результатов:

| Backend | Mode | Size | Workers / ranks x threads | Median time, s | Speedup vs seq | Efficiency | Notes |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| SEQ | task_run | `15000 x 15000`, density `0.001` | `1` | 0.1341834602 | 1.000 | 1.000 | baseline |
| OMP | task_run | `15000 x 15000`, density `0.001` | `1` | 0.2117667000 | 0.634 | 0.634 | overhead выше baseline |
| OMP | task_run | `15000 x 15000`, density `0.001` | `2` | 0.1487973600 | 0.902 | 0.451 | умеренный выигрыш |
| OMP | task_run | `15000 x 15000`, density `0.001` | `4` | 0.1158729000 | 1.158 | 0.290 | лучший диапазон для OMP |
| OMP | task_run | `15000 x 15000`, density `0.001` | `8` | 0.1188691800 | 1.129 | 0.141 | насыщение |
| OMP | task_run | `15000 x 15000`, density `0.001` | `12` | 0.1311886200 | 1.023 | 0.085 | выигрыш почти исчез |
| OMP | task_run | `15000 x 15000`, density `0.001` | `16` | 0.1457355600 | 0.921 | 0.058 | деградация |
| TBB | task_run | `15000 x 15000`, density `0.001` | `1` | 0.1847540202 | 0.726 | 0.726 | близко к baseline |
| TBB | task_run | `15000 x 15000`, density `0.001` | `2` | 0.1124937002 | 1.193 | 0.596 | хороший баланс |
| TBB | task_run | `15000 x 15000`, density `0.001` | `4` | 0.0786577202 | 1.706 | 0.426 | сильный рост |
| TBB | task_run | `15000 x 15000`, density `0.001` | `8` | 0.0699605602 | 1.918 | 0.240 | лучший диапазон |
| TBB | task_run | `15000 x 15000`, density `0.001` | `12` | 0.0697641002 | 1.923 | 0.160 | насыщение |
| TBB | task_run | `15000 x 15000`, density `0.001` | `16` | 0.0698530402 | 1.921 | 0.120 | дополнительный выигрыш минимален |
| STL | task_run | `15000 x 15000`, density `0.001` | `1` | 0.2249369402 | 0.597 | 0.597 | ручной runtime проигрывает baseline |
| STL | task_run | `15000 x 15000`, density `0.001` | `2` | 0.1386025002 | 0.968 | 0.484 | почти выход на baseline |
| STL | task_run | `15000 x 15000`, density `0.001` | `4` | 0.1100570602 | 1.219 | 0.305 | заметный выигрыш |
| STL | task_run | `15000 x 15000`, density `0.001` | `8` | 0.0987977002 | 1.358 | 0.170 | лучший диапазон STL |
| STL | task_run | `15000 x 15000`, density `0.001` | `12` | 0.1027455802 | 1.306 | 0.109 | overhead уже заметен |
| STL | task_run | `15000 x 15000`, density `0.001` | `16` | 0.1199071202 | 1.119 | 0.070 | create/join overhead |
| ALL | task_run | `15000 x 15000`, density `0.001` | `1 x 1` | 0.3087566400 | 0.435 | 0.435 | hybrid overhead без выигрыша |
| ALL | task_run | `15000 x 15000`, density `0.001` | `2 x 2` | 0.1818735600 | 0.738 | 0.184 | коммуникации уже ощутимы |
| ALL | task_run | `15000 x 15000`, density `0.001` | `4 x 4` | 0.1595650400 | 0.841 | 0.053 | лучший медианный режим, но хуже SEQ |
| ALL | task_run | `15000 x 15000`, density `0.001` | `8 x 8` | 0.2955727200 | 0.454 | 0.007 | communication-bound |
| ALL | task_run | `15000 x 15000`, density `0.001` | `12 x 12` | 0.4771294800 | 0.281 | 0.002 | сильная деградация |
| ALL | task_run | `15000 x 15000`, density `0.001` | `16 x 16` | 1.0694690600 | 0.125 | 0.000 | сборка и обмены доминируют |

## 6. Интерпретация различий

`SEQ` показывает реальную цену базового алгоритма без runtime-накладных расходов,
поэтому именно его время используется знаменателем в формулах ускорения.

`OMP` дает только умеренное ускорение. Лучший медианный результат достигается
около `4` потоков, а дальше нерегулярный доступ к строкам `B`, memory-bound
характер нагрузки и последовательная финализация CRS начинают перекрывать пользу
от роста числа потоков.

`TBB` оказался самым сильным backend-ом на этой задаче. Сочетание
`parallel_for`, `blocked_range`, удачного `grain_size` и thread-local workspace
дает наилучший выигрыш, но после диапазона `8–12` workers ускорение уже
практически насыщается.

`STL` дает положительное ускорение, но уступает `TBB`. Лучший результат
достигается около `8` потоков; дальше сильнее проявляются накладные расходы
`create/join` и цена ручного управления балансировкой.

`ALL` на этой машине и на этом размере задачи не показывает реального ускорения
относительно `SEQ`. Даже лучшая медианная конфигурация `4 x 4` остается медленнее
baseline. Причина в том, что к локальному вычислению добавляются `Broadcast`,
`Scatterv`, `Gatherv`, упаковка комплексных значений и финальная сборка на
`rank 0`, из-за чего схема быстро становится communication-bound.

## 7. Репродуцируемость

Команды сборки:

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build -D USE_FUNC_TESTS=ON -D USE_PERF_TESTS=ON -D CMAKE_BUILD_TYPE=Release
cmake -S . -B build -D ENABLE_ADDRESS_SANITIZER=ON -D CMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --config Release --parallel
```

Команды запуска functional-тестов:

```powershell
cd build\bin
.\ppc_func_tests.exe --gtest_filter="*ermakov_a_spar_mat_mult_seq_enabled*"
$env:OMP_NUM_THREADS='8'; .\ppc_func_tests.exe --gtest_filter="*ermakov_a_spar_mat_mult_omp_enabled*"
$env:PPC_NUM_THREADS='8'; .\ppc_func_tests.exe --gtest_filter="*ermakov_a_spar_mat_mult_tbb_enabled*"
$env:PPC_NUM_THREADS='8'; .\ppc_func_tests.exe --gtest_filter="*ermakov_a_spar_mat_mult_stl_enabled*"
$env:PPC_NUM_THREADS='8'; mpiexec -n 8 .\ppc_func_tests.exe --gtest_filter="*ermakov_a_spar_mat_mult_all_enabled*"
```

Команды получения основных `task_run`-замеров:

```powershell
.\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_seq_enabled"
$env:OMP_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_omp_enabled"
$env:PPC_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_tbb_enabled"
$env:PPC_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_stl_enabled"
$env:PPC_NUM_THREADS='8'; mpiexec -n 8 .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_all_enabled"
```

## 8. Заключение

Все пять реализаций корректно решают задачу умножения `CRS x CRS -> CRS` и
проходят общий набор functional-тестов. По медианам `1 прогрев + 5 запусков`
лучшей по ускорению оказалась `TBB`; `STL` дал положительное, но более слабое
ускорение; `OMP` показал только умеренный выигрыш; гибридная `ALL` методически
корректна, но на данной машине не опережает baseline `SEQ`.

Ограничения сравнения:

- измерения выполнены на одной машине и в одном окружении;
- в отчете приведены медианы по `5` измеренным запускам после `1` прогревочного;
- для `ALL` эффективность нормировалась по
  `total_workers = ranks * threads_per_rank`.

Что можно улучшить дальше:

- уменьшить стоимость финальной сборки CRS;
- снизить цену коммуникаций в `ALL`;
- глубже поработать с локальностью данных и балансировкой нагрузки.

## 9. Источники

1. Материалы и инфраструктура курса Parallel Programming Course.
2. OpenMP Architecture Review Board. OpenMP API Specification.
3. MPI Forum. MPI: A Message-Passing Interface Standard.
4. UXL Foundation. oneTBB documentation.
5. cppreference.com: `std::thread`, `std::vector`, `std::complex`.
6. Microsoft MPI documentation.

## 10. Приложение

Локальные технические детали вынесены в подотчеты:

- `seq/report.md`
- `omp/report.md`
- `tbb/report.md`
- `stl/report.md`
- `all/report.md`
