# Умножение плотных матриц. Элементы типа `double`. Алгоритм Штрассена

- Student: Мухаммадхон Исрам
- Variant: 3
- Local reports: `seq/report.md`, `omp/report.md`, `tbb/report.md`, `stl/report.md`, `all/report.md`

## 1. Введение

В работе реализованы и сравниваются пять версий умножения плотных квадратных матриц типа `double`
алгоритмом Штрассена: последовательная (`SEQ`), с OpenMP (`OMP`), с Intel TBB (`TBB`),
со стандартной библиотекой C++ (`STL`) и гибридная (`ALL`, MPI + OMP).

Алгоритм Штрассена хорошо подходит для сравнения моделей параллелизма: на верхнем уровне
есть семь независимых произведений `M1..M7`, которые удобно запускать параллельно,
а базовое блочное умножение при малых размерах даёт предсказуемую нагрузку.
Версия `ALL` использует гибридную схему `MPI + OMP` и рассматривается отдельно
от чисто потоковых backend-ов.

## 2. Единая постановка задачи

Требуется вычислить:

```text
C = A * B
```

где `A` и `B` — квадратные матрицы размера `n x n` с элементами `double`.

Вход: размер `n`, матрица `A` (`n*n` элементов), матрица `B` (`n*n` элементов).

Выход: матрица `C` размера `n x n`.

Ограничения:

- `n > 0`, размеры входных векторов строго равны `n*n`;
- если `n` не степень двойки, матрицы дополняются нулями до ближайшей степени двойки;
- после вычисления padding удаляется;
- критерий корректности: совпадение с эталонным обычным умножением с допуском `1e-9`.

## 3. Единая методика эксперимента

Локальная среда:

- CPU: Intel Core i5-10300H, 4 ядра / 8 потоков;
- RAM: 8 GB;
- OS: Fedora Linux (виртуальная машина);
- compiler: GCC, build type: `Release`.

Команды сборки:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Функциональные тесты:

```bash
cd build/bin
mpiexec -n 1 ./ppc_func_tests --gtest_filter=*muhammadkhon_i_stressen_alg*
```

Performance-тесты:

```bash
cd build/bin
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg*
```

Что важно по методике:

- функциональные тесты и performance в отчёте не смешиваются;
- в сравнительную таблицу вынесены значения `task_run`;
- режим `pipeline` тоже запускался, но в отдельную сравнительную строку не включался;
- текущие числа — это контрольные локальные замеры, без серии повторов с медианой.

## 4. Сводка корректности

Для всех backend-ов функциональные тесты проверяют результат относительно
эталонного обычного умножения матриц с допуском `1e-9`.

Проверялись размеры:

- `3x3`, `4x4` — простейшие случаи и степень двойки;
- `15x5x15`, `16x16` — прямоугольные и средние размеры;
- `63x63`, `64x64` — граница kCutoff и степень двойки.

Локально:

- `SEQ`, `OMP`, `TBB`, `STL` — все функциональные тесты пройдены;
- `ALL` — пройдены на `mpiexec -n 2`.

## 5. Агрегированные результаты

| Backend | Workers / config | Time, s | Speedup vs SEQ | Efficiency | Notes |
| --- | --- | ---: | ---: | ---: | --- |
| seq | 1 | 1.124371 | 1.00 | N/A | baseline |
| omp | 8 threads | 0.441205 | 2.55 | 31.9% | thread-based |
| tbb | 8 workers | 0.398617 | 2.82 | 35.2% | thread-based |
| stl | 8 workers | 0.412883 | 2.72 | 34.0% | thread-based |
| all | 1 rank x 8 threads | 0.096114 | 11.70 | 146.2% | hybrid `MPI + OMP`, `mpiexec -n 1` |

Формулы:

```text
Speedup = T_seq / T_backend
Efficiency = Speedup / Workers * 100%
```

## 6. Интерпретация различий

### SEQ

`SEQ` — это baseline. Его задача не в скорости, а в том, чтобы быть
эталоном по корректности и отправной точкой для расчёта ускорения.
Алгоритм реализован через `std::function impl` внутри обёртки `StrassenSeq`,
что позволило избежать рекурсии свободных функций и предупреждений clang-tidy.

### OMP

`OMP` дал хороший результат среди потоковых backend-ов. Верхний уровень
`M1..M7` параллелится через `#pragma omp task` внутри `#pragma omp single`,
а базовое блочное умножение работает последовательно.
`default(none)` заставил явно перечислять все переменные, что
сделало код понятнее. При одном потоке версия ведёт себя так же, как baseline.

### TBB

`TBB` показал лучшее время среди потоковых backend-ов. Верхний уровень
`M1..M7` запускается через `oneapi::tbb::parallel_invoke` — каждая ветвь
оформлена как отдельная лямбда. `parallel_for` и `blocked_range`
не использовались, потому что задачи дискретные.
Число потоков ограничивается через `global_control::max_allowed_parallelism`.

### STL

`STL` тоже дал ускорение, близкое к OMP. Верхний уровень `M1..M7`
запускается через `std::async(std::launch::async, ...)`, синхронизация — через
`future.get()`. Каждая задача пишет в свой буфер, гонок нет, mutex не нужен.
Версия оказалась чуть медленнее TBB, что ожидаемо: нет автоматической
балансировки, только семь фиксированных задач.

### ALL

`ALL` — гибридная версия `MPI + OMP`. Семь подзадач `M1..M7` распределяются
между MPI-процессами по схеме round-robin, каждый процесс считает свои задачи
через OMP, вклады собираются через `MPI_Reduce`. Рассматривать эту версию
нужно отдельно: важна конфигурация `ranks x threads`, а не просто
количество воркеров.

## 7. Репродуцируемость

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

cd build/bin

# Функциональные тесты
mpiexec -n 1 ./ppc_func_tests --gtest_filter=*muhammadkhon_i_stressen_alg*
mpiexec -n 2 ./ppc_func_tests --gtest_filter=*muhammadkhon_i_stressen_alg_all*

# Performance
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg_seq*
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg_omp*
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg_tbb*
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg_stl*
mpiexec -n 1 ./ppc_perf_tests --gtest_filter=*muhammadkhon_i_stressen_alg_all*
```

Практические ограничения:

- замеры выполнялись на виртуальной машине, что может влиять на абсолютные числа;
- без серии повторов с медианой результаты стоит рассматривать как ориентировочные;
- для `ALL` конфигурация `ranks x threads` важна отдельно.

## 8. Заключение

Из всех потоковых backend-ов лучшее локальное время показал `TBB`.
Это ожидаемо: `parallel_invoke` хорошо подходит для семи независимых задач
верхнего уровня Штрассена, а work-stealing планировщик балансирует нагрузку
без ручного вмешательства.

`OMP` получился чуть медленнее, но зато проще в реализации и понятнее
по структуре директив. `STL` занял промежуточную позицию между ними.

`ALL` — отдельный гибридный backend, его нужно рассматривать
со своей схемой измерений для разных `ranks x threads`.

Итог:

- `SEQ` — честный baseline;
- `OMP`, `TBB`, `STL` — понятное локальное ускорение ~2.5–2.8x на 8 потоках;
- `ALL` — гибридный backend, корректно работает на `mpiexec -n 2`.

## 9. Источники

1. Документация курса Parallel Programming Course, ННГУ.
2. V. Strassen, *Gaussian elimination is not optimal*, Numerische Mathematik, 1969.
3. [OpenMP Specifications](https://www.openmp.org/specifications/)
4. [MPI Forum Documentation](https://www.mpi-forum.org/docs/)
5. [oneTBB Documentation](https://uxlfoundation.github.io/oneTBB/)

## 10. Приложение

Ключевые файлы:

- `tasks/muhammadkhon_i_stressen_alg/seq/src/ops_seq.cpp`
- `tasks/muhammadkhon_i_stressen_alg/omp/src/ops_omp.cpp`
- `tasks/muhammadkhon_i_stressen_alg/tbb/src/ops_tbb.cpp`
- `tasks/muhammadkhon_i_stressen_alg/stl/src/ops_stl.cpp`
- `tasks/muhammadkhon_i_stressen_alg/all/src/ops_all.cpp`
- `tasks/muhammadkhon_i_stressen_alg/tests/functional/main.cpp`
- `tasks/muhammadkhon_i_stressen_alg/tests/performance/main.cpp`
