# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Student: Зорин Данила Артёмович
- Variant: 3
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`,
  `stl/report.md`, `all/report.md`

## 1. Введение

В работе сравниваются пять реализаций одной и той же задачи:
`SEQ`, `OMP`, `TBB`, `STL` и `ALL`.

Во всех вариантах вычисляется произведение плотных квадратных матриц
типа `double` алгоритмом Штрассена. Эта задача удобна для сравнения разных
моделей параллелизма, потому что у неё есть независимые ветви `M1..M7`
и заметные накладные расходы на временные буферы, padding и сборку результата.

В версии `ALL` используется гибридная схема `MPI + OMP`, поэтому в сводном
отчёте её корректнее рассматривать отдельно от чисто потоковых backend-ов.

## 2. Единая постановка задачи

Требуется вычислить:

```text
C = A * B
```

где `A` и `B` — квадратные матрицы размера `n x n` с элементами `double`.

Вход:

- размер `n`;
- матрица `A`;
- матрица `B`.

Выход:

- матрица `C`.

Ограничения:

- `n > 0`;
- размеры `A` и `B` должны совпадать с `n * n`;
- если `n` не степень двойки, матрицы дополняются нулями до ближайшей
  степени двойки;
- после вычисления padding удаляется.

Критерий корректности:

- результат сверяется с обычным эталонным умножением матриц;
- сравнение выполняется с допуском `1e-9`.

## 3. Единая методика эксперимента

Локальная среда:

- CPU: AMD Ryzen 5 2600, 6 ядер / 12 потоков;
- RAM: 16 GB;
- OS: Windows 11 x64;
- compiler: MSVC / Visual Studio 2022;
- build type: `Release`.

Команды сборки:

```powershell
cmake -S . -B build_release
cmake --build build_release --config Release --parallel
```

Функциональные тесты:

```powershell
cd build\bin
mpiexec -n 1 .\ppc_func_tests --gtest_filter=*zorin_d_strassen_alg_matrix_seq*
```

Performance-тесты:

```powershell
cd build\bin
mpiexec -n 1 .\ppc_perf_tests --gtest_filter=StrassenPerf/ZorinDRunPerfTests.*
```

Что важно по методике:

- функциональные тесты и performance в отчёте не смешиваются;
- в сравнительную таблицу вынесены значения `task_run`;
- режим `pipeline` тоже запускался, но в отдельную сравнительную строку
  я его не смешивал с `task_run`;
- текущие числа — это контрольные локальные замеры, а не полноценная серия
  повторов с медианой.

## 4. Сводка корректности

В `tests/functional/main.cpp` для всех backend-ов генерируются одинаковые
входные матрицы, после чего результат сверяется с эталонным обычным
умножением.

Проверялись размеры:

- `1`;
- `2`;
- `3`;
- `5`;
- `8`;
- `9`;
- `16`.

Эти размеры нужны не случайно:

- `1` и `2` проверяют простейшие случаи;
- `3`, `5` и `9` проверяют корректность padding;
- `8` и `16` дают уже типичные размеры для рекурсивного алгоритма.

Локально у меня:

- `ALL` проходит на `mpiexec -n 2`;
- полный набор `zorin_d_strassen_alg_matrix_seq` на `mpiexec -n 1`
  проходит полностью.

## 5. Агрегированные результаты

| Backend | Workers / config | Time, s | Speedup vs SEQ | Efficiency | Notes |
| --- | --- | ---: | ---: | ---: | --- |
| seq | 1 | 1.023319 | 1.00 | N/A | baseline |
| omp | 8 threads | 0.420932 | 2.43 | 30.4% | thread-based |
| tbb | 8 workers | 0.482507 | 2.12 | 26.5% | thread-based |
| stl | 8 workers | 0.491042 | 2.08 | 26.0% | thread-based |
| all | 1 rank x 8 threads | 0.088702 | 11.54 | 144.3% | hybrid `MPI + OMP`, measured at `mpiexec -n 1` |

Формулы:

```text
Speedup = T_seq / T_backend
Efficiency = Speedup / Workers * 100%
```


## 6. Интерпретация различий

### SEQ

`SEQ` — это baseline. Его задача не в скорости, а в том, чтобы быть
эталоном по корректности и отправной точкой для расчёта ускорения.

### OMP

`OMP` дал лучший результат среди чисто потоковых backend-ов.
Для этой задачи это выглядит логично: верхние ветви Штрассена хорошо
ложатся на `omp task`, а код остаётся достаточно простым.

### TBB

`TBB` показал стабильное ускорение и оказался близким к `OMP`,
но всё же немного медленнее на моих локальных замерах.

### STL

`STL` тоже ускорился заметно. По времени он получился близок к `TBB`,
хотя я ожидал, что из-за более ручного управления задачами он проиграет сильнее.

### ALL

`ALL` здесь нужно рассматривать как отдельный гибридный backend.

## 7. Репродуцируемость

Минимальный набор команд для воспроизведения:

```powershell
cmake -S . -B build_release
cmake --build build_release --config Release --parallel

cd build\bin
mpiexec -n 1 .\ppc_func_tests --gtest_filter=*zorin_d_strassen_alg_matrix_seq*
mpiexec -n 2 .\ppc_func_tests --gtest_filter=*zorin_d_strassen_alg_matrix_seq_all_enabled*
mpiexec -n 1 .\ppc_perf_tests --gtest_filter=StrassenPerf/ZorinDRunPerfTests.*
```

Практическое ограничение:

- локальные замеры выполнялись на Windows;
- поведение в CI может отличаться, особенно на `macOS`;
- без отдельной серии повторов и замеров для разных `ranks x threads`
  итоговые числа для `ALL` стоит считать предварительными.

## 8. Заключение

Если смотреть только на thread-based реализации, самым удачным вариантом у меня
получился `OMP`. Он дал лучшее локальное время и при этом остался довольно
понятным по структуре.

В версии `ALL` используется гибридная схема `MPI + OMP`, и функциональные
тесты для неё проходят корректно и на `mpiexec -n 2`.

Итог по работе для меня такой:

- `SEQ` нужен как честный baseline;
- `OMP`, `TBB` и `STL` дают понятное локальное ускорение;
- `ALL` стоит оценивать как отдельный гибридный backend
  со своей схемой измерений.

## 9. Источники

1. Документация курса Parallel Programming Course.
2. V. Strassen, *Gaussian elimination is not optimal*,
   Numerische Mathematik, 1969.
3. [OpenMP Specifications](https://www.openmp.org/specifications/)
4. [MPI Forum Documentation](https://www.mpi-forum.org/docs/)
5. [oneTBB Documentation](https://uxlfoundation.github.io/oneTBB/)

## 10. Приложение

Ключевые файлы:

- `tasks/zorin_d_strassen_alg_matrix_seq/seq/src/ops_seq.cpp`
- `tasks/zorin_d_strassen_alg_matrix_seq/omp/src/ops_omp.cpp`
- `tasks/zorin_d_strassen_alg_matrix_seq/tbb/src/ops_tbb.cpp`
- `tasks/zorin_d_strassen_alg_matrix_seq/stl/src/ops_stl.cpp`
- `tasks/zorin_d_strassen_alg_matrix_seq/all/src/ops_all.cpp`
- `tasks/zorin_d_strassen_alg_matrix_seq/tests/functional/main.cpp`
- `tasks/zorin_d_strassen_alg_matrix_seq/tests/performance/main.cpp`
