# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – строковый (CRS)

- Студент: Ермаков Алексей Викторович
- Группа 3823Б1ПР3
- Вариант: 6
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

В работе рассматривается задача умножения двух разреженных комплексных
матриц в формате CRS. Эта задача удобна для сравнения разных моделей
параллелизма, поскольку вычислительное ядро у всех реализаций общее, а
различаются способы организации исполнения: последовательный baseline,
OpenMP, oneTBB, `std::thread` и гибридная схема `MPI + OMP`.

## 2. Единая постановка задачи

Вход:

- две матрицы `A` и `B` в формате CRS;
- каждая матрица описывается полями `rows`, `cols`, `values`, `col_index`, `row_ptr`.

Выход:

- матрица `C = A * B` в формате CRS.

Критерии корректности:

- `A.cols == B.rows`;
- `row_ptr.size() == rows + 1`;
- `values.size() == col_index.size()`;
- `row_ptr.front() == 0`;
- `row_ptr.back() == nnz`;
- все индексы `col_index[k]` лежат в диапазоне `[0, cols)`;
- результат совпадает с плотным эталоном умножения.

## 3. Единая методика эксперимента

Окружение:

- CPU: AMD Ryzen 5 1600
- RAM: 8 GB
- OS: Windows 11
- compiler: MSVC
- build type: `Release`
- MPI: MS-MPI

Входные данные:

- functional-тесты: `3x3`, `10x10`, `20x20`, `30x30`;
- performance-тест: `15000 x 15000`, плотность `0.001`;
- генерация performance-входа детерминированная, через фиксированные seed-ы `0x13579BDFU` и `0x2468ACE0U`.

Методика замера:

- сравнение ведется по `task_run`;
- для каждой конфигурации выполняется `1` прогревочный запуск;
- затем выполняется `5` измеренных запусков;
- в таблицы выносится медиана этих `5` запусков.

Формулы:

- `speedup = T_seq / T_backend`;
- для `OMP`, `TBB`, `STL`: `efficiency = speedup / threads`;
- для `ALL`: `efficiency = speedup / total_workers`, где `total_workers = ranks * threads_per_rank`.

## 4. Сводка корректности

Все backend-ы проверялись общей тестовой инфраструктурой проекта:

- functional-набор из `tests/functional/main.cpp`;
- performance-набор из `tests/performance/main.cpp`.

Корректность подтверждалась:

- сравнением с baseline `SEQ`;
- сравнением с плотным эталоном `DenseMul`;
- проверкой CRS-структуры результата.

Для `ALL` в functional-режиме итоговая матрица после сборки на `rank 0`
дополнительно рассылается всем rank-ам, поэтому `PostProcessingImpl()`
получает согласованный результат на каждом процессе.

## 5. Агрегированные результаты

Последовательный baseline:

- `T_seq = 0.1300215802 c`.

Сводная таблица `task_run`:

| Backend | Mode | Size | Workers / ranks x threads | Median time, s | Speedup vs seq | Efficiency | Notes |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| SEQ | task_run | `15000 x 15000`, density `0.001` | `1` | 0.1300215802 | 1.000 | 1.000 | baseline |
| OMP | task_run | `15000 x 15000`, density `0.001` | `1` | 0.2117667000 | 0.614 | 0.614 | overhead выше baseline |
| OMP | task_run | `15000 x 15000`, density `0.001` | `2` | 0.1487973600 | 0.874 | 0.437 | умеренный выигрыш |
| OMP | task_run | `15000 x 15000`, density `0.001` | `4` | 0.1158729000 | 1.122 | 0.280 | лучший диапазон для OMP |
| OMP | task_run | `15000 x 15000`, density `0.001` | `8` | 0.1188691800 | 1.094 | 0.137 | насыщение |
| OMP | task_run | `15000 x 15000`, density `0.001` | `12` | 0.1311886200 | 0.991 | 0.083 | выигрыш почти исчез |
| OMP | task_run | `15000 x 15000`, density `0.001` | `16` | 0.1457355600 | 0.892 | 0.056 | деградация |
| TBB | task_run | `15000 x 15000`, density `0.001` | `1` | 0.1847540202 | 0.704 | 0.704 | близко к baseline |
| TBB | task_run | `15000 x 15000`, density `0.001` | `2` | 0.1124937002 | 1.156 | 0.578 | хороший баланс |
| TBB | task_run | `15000 x 15000`, density `0.001` | `4` | 0.0786577202 | 1.653 | 0.413 | сильный рост |
| TBB | task_run | `15000 x 15000`, density `0.001` | `8` | 0.0699605602 | 1.859 | 0.232 | лучший диапазон |
| TBB | task_run | `15000 x 15000`, density `0.001` | `12` | 0.0697641002 | 1.864 | 0.155 | насыщение |
| TBB | task_run | `15000 x 15000`, density `0.001` | `16` | 0.0698530402 | 1.861 | 0.116 | дополнительный выигрыш минимален |
| STL | task_run | `15000 x 15000`, density `0.001` | `1` | 0.2249369402 | 0.578 | 0.578 | runtime дороже baseline |
| STL | task_run | `15000 x 15000`, density `0.001` | `2` | 0.1386025002 | 0.938 | 0.469 | почти уровень baseline |
| STL | task_run | `15000 x 15000`, density `0.001` | `4` | 0.1100570602 | 1.181 | 0.295 | заметный выигрыш |
| STL | task_run | `15000 x 15000`, density `0.001` | `8` | 0.0987977002 | 1.316 | 0.164 | лучший диапазон STL |
| STL | task_run | `15000 x 15000`, density `0.001` | `12` | 0.1027455802 | 1.266 | 0.106 | overhead уже заметен |
| STL | task_run | `15000 x 15000`, density `0.001` | `16` | 0.1199071202 | 1.084 | 0.068 | create/join overhead |
| ALL | task_run | `15000 x 15000`, density `0.001` | `1 x 1` | 0.1952643600 | 0.666 | 0.666 | MPI-обвязка без выигрыша |
| ALL | task_run | `15000 x 15000`, density `0.001` | `2 x 2` | 0.1071202600 | 1.214 | 0.303 | быстрее SEQ |
| ALL | task_run | `15000 x 15000`, density `0.001` | `4 x 4` | 0.0903336400 | 1.439 | 0.090 | лучший режим ALL |
| ALL | task_run | `15000 x 15000`, density `0.001` | `8 x 8` | 0.1192515800 | 1.090 | 0.017 | еще быстрее SEQ |
| ALL | task_run | `15000 x 15000`, density `0.001` | `12 x 12` | 0.1414913600 | 0.919 | 0.006 | коммуникации доминируют |
| ALL | task_run | `15000 x 15000`, density `0.001` | `16 x 16` | 0.2347756400 | 0.554 | 0.002 | сильная деградация |

## 6. Интерпретация различий

`SEQ` задает baseline для сравнения по времени и корректности.

`OMP` дает умеренный выигрыш, но быстро упирается в memory-bound характер
вычислений и последовательную фазу формирования итогового CRS.

`TBB` показывает лучший общий результат. Сочетание `parallel_for`,
`blocked_range` и thread-local буферов дает наиболее стабильное ускорение
на этой задаче.

`STL` тоже ускоряет вычисление, но уступает `TBB` из-за более высокой цены ручного управления потоками.

`ALL` реализует двухуровневую схему:

- строки `A` распределяются между MPI-процессами;
- полная матрица `B` рассылается всем rank-ам;
- локальный блок строк внутри rank-а обрабатывается через OpenMP;
- глобальный CRS собирается на `rank 0`.

В этой реализации:

- границы блоков строк строятся по оценке стоимости умножения `BuildRowBounds(a, b, ...)`;
- комплексные значения передаются через собственный MPI-тип напрямую;
- число OpenMP-потоков внутри rank-а ограничивается с учетом `hardware_concurrency()`.

На практике `ALL` дает выигрыш относительно `SEQ` в конфигурациях `2 x 2`,
`4 x 4` и `8 x 8`, а лучший результат достигается в `4 x 4`. При
дальнейшем росте числа rank-ов и потоков начинают доминировать стоимость
MPI-коммуникаций, синхронизации и финальной сборки результата на `rank 0`.

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

Команды получения `task_run`:

```powershell
.\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_seq_enabled"
$env:OMP_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_omp_enabled"
$env:PPC_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_tbb_enabled"
$env:PPC_NUM_THREADS='8'; .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_stl_enabled"
$env:PPC_NUM_THREADS='8'; mpiexec -n 8 .\ppc_perf_tests.exe --gtest_filter="*task_run_ermakov_a_spar_mat_mult_all_enabled"
```

## 8. Заключение

Все пять реализаций корректно решают задачу умножения `CRS x CRS -> CRS`.
По производительности лучшим backend-ом в проведенных измерениях остается
`TBB`.

Гибридная версия `ALL` показывает выигрыш относительно последовательного
baseline в умеренных конфигурациях `2 x 2`, `4 x 4` и `8 x 8`. Лучший
режим - `4 x 4`, где достигается наибольший speedup для этой реализации.

Ограничения сравнения:

- замеры выполнены на одной машине и в одном окружении;
- в таблицах приведены медианы `1 прогрев + 5 измерений`;
- для `ALL` эффективность нормировалась по `total_workers = ranks * threads_per_rank`.

Что можно улучшить дальше:

- уменьшить стоимость финальной сборки CRS на `rank 0`;
- сократить цену рассылки полной матрицы `B`;
- глубже подбирать баланс `ranks x threads` под конкретную машину.

## 9. Источники

1. Материалы и инфраструктура курса Parallel Programming Course.
2. OpenMP Architecture Review Board. OpenMP API Specification.
3. MPI Forum. MPI: A Message-Passing Interface Standard.
4. UXL Foundation. oneTBB documentation.
5. cppreference.com: `std::thread`, `std::vector`, `std::complex`.
6. Microsoft MPI documentation.

## 10. Приложение

Подробные технологические детали вынесены в локальные отчеты:

- `seq/report.md`
- `omp/report.md`
- `tbb/report.md`
- `stl/report.md`
- `all/report.md`
