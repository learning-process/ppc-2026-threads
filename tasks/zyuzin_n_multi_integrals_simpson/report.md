# Вычисление многомерных интегралов с использованием многошаговой схемы (метод Симпсона)

- Student: Зюзин Никита Михайлович, группа 3823Б1ПР2
- Variant: 11

Local reports: `seq/report.md`, `omp/report.md`, `stl/report.md`, `tbb/report.md`, `all/report.md`

## 1. Введение

Численное интегрирование многомерных функций является фундаментальной задачей вычислительной
математики. Оно применяется в физическом моделировании, теории вероятностей, машинном обучении
и задачах оптимизации, где аналитическое вычисление интеграла недоступно или нецелесообразно.

Метод Симпсона (составная формула Симпсона) даёт точность порядка `O(h⁴)` по каждой оси при
достаточной гладкости подынтегральной функции. Для d-мерного пространства на основе декартового
произведения одномерных сеток суммарное число вычислений функции равно `∏(n_steps[dim]+1)`, что
при средних `n_steps` ~100 и d=3 даёт ~10⁶ вычислений — объём, хорошо поддающийся параллелизму.

Данная задача особенно удобна для сравнения моделей параллелизма, поскольку:

- Все вычисления на узлах решётки полностью независимы (нет зависимостей по данным);
- Нагрузка на каждый узел равномерна (одинаковое число операций);
- Единственный результат — скалярная сумма — легко агрегируется редукцией.

В работе реализованы и сравнены пять вариантов:

- SEQ — последовательная эталонная версия
- OMP — параллелизм через OpenMP
- STL — ручное управление потоками через `std::thread`
- TBB — параллелизм через oneTBB `parallel_reduce`
- ALL — гибридная версия (MPI + OpenMP)

## 2. Единая постановка задачи

### Входные данные

Структура `SimpsonInput` (`common/include/common.hpp`):

- `lower_bounds` — вектор нижних пределов интегрирования по каждому измерению;
- `upper_bounds` — вектор верхних пределов интегрирования по каждому измерению;
- `n_steps` — вектор числа шагов по каждому измерению;
- `func` — подынтегральная функция `std::function<double(const std::vector<double>&)>`.

### Выходные данные

Приближённое значение интеграла типа `double`.

### Ограничения

- Размеры `lower_bounds`, `upper_bounds`, `n_steps` совпадают и ненулевые;
- `lower_bounds[i] <= upper_bounds[i]` для каждого i;
- `n_steps[i] > 0` и `n_steps[i] % 2 == 0` (требование формулы Симпсона);
- `func` задана (непустой `std::function`).

### Алгоритм (общий для всех реализаций)

1. Вычисление шагов: `h[d] = (upper[d] - lower[d]) / n_steps[d]`
2. Общее число узлов: `total_points = ∏(n_steps[d]+1)`
3. Для каждого узла `point_idx ∈ [0, total_points)`:
   - Декодирование в многомерный индекс: `indices[d] = (point_idx / ∏_{k<d}(n_steps[k]+1)) % (n_steps[d]+1)`
   - Координата: `x[d] = lower[d] + indices[d] * h[d]`
   - Вес: `w = ∏ GetSimpsonWeight(indices[d], n_steps[d])`
   - `sum += w * func(x)`
4. Результат: `sum * ∏(h[d]/3)`

Вес одномерного узла Симпсона:

```text
w(0, n) = 1,  w(n, n) = 1
w(i, n) = 4   если i нечётное
w(i, n) = 2   если i чётное, 0 < i < n
```

### Критерий корректности

`|result - expected| < 1e-3`

## 3. Единая методика эксперимента

### Hardware / OS

- Процессор: AMD Ryzen 5 5600X
- Ядра/потоки: 6 ядер / 12 потоков
- ОЗУ: 16 GB
- ОС: Ubuntu 25.10
- Архитектура: x64

### Toolchain

- Компилятор: GCC 15.2.0
- IDE: Visual Studio Code 2026
- Тип сборки: Release
- Система сборки: CMake
- Версия MPI: Open MPI 5.0.8

### Переменные окружения

| Технология | Переменные                          |
|------------|-------------------------------------|
| SEQ        | (не требуются)                      |
| OMP        | `OMP_NUM_THREADS=N`                 |
| STL        | `PPC_NUM_THREADS=N`                 |
| TBB        | `PPC_NUM_THREADS=N`                 |
| ALL        | `OMP_NUM_THREADS=N`, `mpirun -np P` |

### Входные данные для performance-теста

3D интеграл: `sin(x)·cos(y)·exp(z)` на `[0,π]×[0,π/2]×[0,1]`,
`n_steps = {180, 180, 120}`.

Всего узлов: `181 × 181 × 121 = 3 964 381 ≈ 4M`.
Ожидаемое значение: `2·(e−1) ≈ 3.4366`.

### Метрики

- Ускорение: `Speedup = T_seq / T_parallel`
- Эффективность: `Efficiency = (Speedup / Workers) × 100%`

## 4. Сводка корректности

### Функциональные тесты (12 случаев, допуск `1e-3`)

| ID | Название       | Размерность | Функция               | Область              | Ожидание |
|----|----------------|-------------|-----------------------|----------------------|----------|
| 0  | 1d_linear      | 1D          | x                     | [0,1]                | 0.5      |
| 1  | 1d_quadratic   | 1D          | x^2                   | [0,1]                | 1/3      |
| 2  | 2d_sum         | 2D          | x+y                   | [0,1]^2              | 1.0      |
| 3  | 2d_product     | 2D          | x·y                   | [0,1]^2              | 0.25     |
| 4  | 2d_sum_squares | 2D          | x^2+y^2               | [0,1]^2              | 2/3      |
| 5  | 2d_constant    | 2D          | 1                     | [0,2]×[0,3]          | 6.0      |
| 6  | 3d_constant    | 3D          | 1                     | [0,1]^3              | 1.0      |
| 7  | 3d_sum         | 3D          | x+y+z                 | [0,1]^3              | 1.5      |
| 8  | 1d_sin         | 1D          | sin(x)                | [0,π]                | 2.0      |
| 9  | 2d_sin_cos     | 2D          | sin(x)·cos(y)         | [0,π/2]^2            | 1.0      |
| 10 | 2d_exp         | 2D          | exp(x+y)              | [0,1]^2              | (e−1)^2  |
| 11 | 3d_product     | 3D          | x·y·z                 | [0,1]^3              | 1/8      |

Все реализации (SEQ, OMP, STL, TBB, ALL) проходят все 12 тестов. ALL-тесты
выполняются только под `mpirun` (в одиночном запуске они пропускаются фреймворком с пометкой
`kALL and kMPI tasks are not under mpirun`).

## 5. Различия реализаций

| Реализация | Параллелизм          | Буфер `point`           | Разбиение           | Редукция              |
|------------|----------------------|-------------------------|---------------------|-----------------------|
| SEQ        | нет                  | per-iteration (vector)  | нет                 | накопление в `sum`    |
| OMP        | OpenMP parallel for  | per-iteration (vector)  | статическое (omp)   | `reduction(+: sum)`   |
| STL        | std::thread          | per-thread (vector)     | блочное вручную     | `std::accumulate`     |
| TBB        | parallel_reduce      | per-grain (vector)      | blocked_range(4096) | TBB-редукция          |
| ALL        | MPI + OpenMP         | per-iteration (vector)  | по рангам + static  | `MPI_Allreduce`       |

## 6. Агрегированные результаты

Запуск: `ppc_perf_tests`, Release-сборка, серии прогонов по командам репродуцируемости
(`OMP/STL/TBB`: 2/4/6/8/10/12 потоков, `ALL`: 2/4 процесса).

| Mode       | Workers (MPI×OMP) | Время, с | Ускорение | Эффективность |
|------------|-------------------|----------|-----------|---------------|
| seq (task) | 1×1               | 1.1135   | 1.00      | 100%          |
| seq (pipe) | 1×1               | 1.0734   | 1.00      | 100%          |
| omp (task) | 1×2               | 0.6173   | 1.80      | 90%           |
| omp (pipe) | 1×2               | 0.5562   | 1.93      | 97%           |
| omp (task) | 1×4               | 0.3628   | 3.07      | 77%           |
| omp (pipe) | 1×4               | 0.2847   | 3.77      | 94%           |
| omp (task) | 1×6               | 0.2627   | 4.24      | 71%           |
| omp (pipe) | 1×6               | 0.2047   | 5.24      | 87%           |
| omp (task) | 1×8               | 0.2617   | 4.25      | 53%           |
| omp (pipe) | 1×8               | 0.2114   | 5.08      | 63%           |
| omp (task) | 1×10              | 0.2043   | 5.45      | 54%           |
| omp (pipe) | 1×10              | 0.1946   | 5.52      | 55%           |
| omp (task) | 1×12              | 0.1960   | 5.68      | 47%           |
| omp (pipe) | 1×12              | 0.1937   | 5.54      | 46%           |
| stl (task) | 1×2               | 0.2033   | 5.48      | 274%          |
| stl (pipe) | 1×2               | 0.2297   | 4.67      | 234%          |
| stl (task) | 1×4               | 0.1145   | 9.73      | 243%          |
| stl (pipe) | 1×4               | 0.1062   | 10.11     | 253%          |
| stl (task) | 1×6               | 0.0747   | 14.91     | 248%          |
| stl (pipe) | 1×6               | 0.0781   | 13.75     | 229%          |
| stl (task) | 1×8               | 0.0882   | 12.63     | 158%          |
| stl (pipe) | 1×8               | 0.0856   | 12.54     | 157%          |
| stl (task) | 1×10              | 0.0807   | 13.79     | 138%          |
| stl (pipe) | 1×10              | 0.0784   | 13.69     | 137%          |
| stl (task) | 1×12              | 0.0752   | 14.82     | 124%          |
| stl (pipe) | 1×12              | 0.0753   | 14.25     | 119%          |
| tbb (task) | 1×2               | 0.1972   | 5.65      | 283%          |
| tbb (pipe) | 1×2               | 0.1987   | 5.40      | 270%          |
| tbb (task) | 1×4               | 0.1028   | 10.83     | 271%          |
| tbb (pipe) | 1×4               | 0.1045   | 10.27     | 257%          |
| tbb (task) | 1×6               | 0.0721   | 15.45     | 258%          |
| tbb (pipe) | 1×6               | 0.0738   | 14.54     | 242%          |
| tbb (task) | 1×8               | 0.0728   | 15.30     | 191%          |
| tbb (pipe) | 1×8               | 0.0722   | 14.87     | 186%          |
| tbb (task) | 1×10              | 0.0698   | 15.94     | 159%          |
| tbb (pipe) | 1×10              | 0.0686   | 15.65     | 156%          |
| tbb (task) | 1×12              | 0.0677   | 16.46     | 137%          |
| tbb (pipe) | 1×12              | 0.0667   | 16.10     | 134%          |
| all (task) | 2×1               | 0.3328   | 3.20      | 160%          |
| all (pipe) | 2×1               | 0.3243   | 3.30      | 165%          |
| all (task) | 2×2               | 0.2710   | 3.92      | 98%           |
| all (pipe) | 2×2               | 0.2711   | 3.94      | 99%           |
| all (task) | 2×3               | 0.2748   | 3.87      | 64%           |
| all (pipe) | 2×3               | 0.2712   | 3.94      | 66%           |
| all (task) | 2×4               | 0.2700   | 3.94      | 49%           |
| all (pipe) | 2×4               | 0.2698   | 3.96      | 50%           |
| all (task) | 2×5               | 0.2763   | 3.85      | 39%           |
| all (pipe) | 2×5               | 0.2759   | 3.87      | 39%           |
| all (task) | 2×6               | 0.2700   | 3.94      | 33%           |
| all (pipe) | 2×6               | 0.2696   | 3.96      | 33%           |
| all (task) | 4×1               | 0.1725   | 6.17      | 154%          |
| all (pipe) | 4×1               | 0.1745   | 6.12      | 153%          |
| all (task) | 4×2               | 0.1376   | 7.73      | 97%           |
| all (pipe) | 4×2               | 0.1381   | 7.74      | 97%           |
| all (task) | 4×3               | 0.1384   | 7.69      | 64%           |
| all (pipe) | 4×3               | 0.1381   | 7.74      | 64%           |

## 7. Интерпретация различий

### SEQ

Последовательная версия — эталон корректности и точка отсчёта производительности. Основной
источник накладных расходов — выделение двух `std::vector` на каждую из ~4M итераций.

### OMP

Минимальные изменения кода: добавлена одна прагма. `reduction(+: sum)` корректно агрегирует
частичные суммы. OMP масштабируется до 12 потоков умеренно: прирост после
6 потоков замедляется, что согласуется с накладными расходами на per-iteration аллокацию `vector`.

### STL

Явное блочное разбиение и переиспользование `point`-буфера на поток. Нет зависимости от
OpenMP или TBB — только стандартная библиотека. На новых измерениях STL показывает высокую
производительность, но с немонотонным поведением (локальная просадка около 8 потоков).

### TBB

`parallel_reduce` с `grain_size = 4096` и work-stealing планировщиком. Буфер `point`
переиспользуется в рамках одного grain (4096 итераций), что снижает давление на аллокатор.
На серии 2..12 потоков TBB показывает наиболее устойчивый рост и лучший результат на 12 потоках.

### ALL

Два уровня параллелизма. MPI делит `total_points` между процессами, OMP параллелизует
локальный поддиапазон. `MPI_Allreduce` — симметричная операция без выделенного мастер-ранга.
В диапазоне до 12 workers лучший режим — `4×2` (близко к нему `4×3`), тогда как
режимы `2×N` после `N=2` масштабируются слабее.

## 8. Репродуцируемость

Команды запуска функциональных тестов

```cpp
// SEQ
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*seq*"

// OMP (2 потока)
$env:OMP_NUM_THREADS=2
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*omp*"

// TBB (2 потока)
$env:PPC_NUM_THREADS=2
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*tbb*"

// STL (2 потока)
$env:PPC_NUM_THREADS=2
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*stl*"

// ALL (2 процесса × 2 потока внутри)
$env:OMP_NUM_THREADS=2
mpiexec -np 2 ./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*all*"
```

Команды запуска тестов производительности

```cpp
// SEQ (baseline)
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*seq*"

// OMP с 2, 4, 6, 8, 10 и 12 потоками
$env:OMP_NUM_THREADS=2
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*omp*"

$env:OMP_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*omp*"

$env:OMP_NUM_THREADS=6
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*omp*"

$env:OMP_NUM_THREADS=8
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*omp*"

$env:OMP_NUM_THREADS=10
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*omp*"

$env:OMP_NUM_THREADS=12
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*omp*"

// TBB с 2, 4, 6, 8, 10 и 12 потоками
$env:PPC_NUM_THREADS=2
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*tbb*"

$env:PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*tbb*"

$env:PPC_NUM_THREADS=6
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*tbb*"

$env:PPC_NUM_THREADS=8
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*tbb*"

$env:PPC_NUM_THREADS=10
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*tbb*"

$env:PPC_NUM_THREADS=12
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*tbb*"

// STL с 2, 4, 6, 8, 10 и 12 потоками
$env:PPC_NUM_THREADS=2
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*stl*"

$env:PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*stl*"

$env:PPC_NUM_THREADS=6
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*stl*"

$env:PPC_NUM_THREADS=8
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*stl*"

$env:PPC_NUM_THREADS=10
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*stl*"

$env:PPC_NUM_THREADS=12
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*stl*"

// ALL с 2 и 4 процессами (по 1 потоку внутри)
$env:OMP_NUM_THREADS=1
mpiexec -np 2 ./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*all*"
mpiexec -np 4 ./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*all*"
```

## 9. Заключение

Реализованы и протестированы пять вариантов вычисления многомерных интегралов методом
Симпсона. Все реализации корректны — проходят полный набор функциональных тестов.

Итоги по реализациям:

- **SEQ** — эталон корректности, без параллелизма, per-iteration выделение `vector`;
- **OMP** — минимальное изменение SEQ, одна прагма, статическое планирование, `reduction`;
- **STL** — ручное управление потоками, переиспользование буфера на поток, `std::accumulate`;
- **TBB** — `parallel_reduce` с `grain_size = 4096`, переиспользование буфера на grain;
- **ALL** — MPI + OpenMP, симметричная сборка через `MPI_Allreduce`, два уровня параллелизма.

По результатам измерений:

- OMP масштабируется до 12 потоков, но с заметным замедлением прироста после 6 потоков;
- STL и TBB существенно быстрее OMP на данной задаче благодаря переиспользованию буфера `point`;
- TBB показывает наиболее стабильное масштабирование на диапазоне 2..12 потоков;
- ALL даёт лучший результат в конфигурациях с 4 MPI-процессами (`4×2` и `4×3`), что подтверждает
  эффективность гибридной схемы при правильном балансе между MPI и OMP уровнями.
