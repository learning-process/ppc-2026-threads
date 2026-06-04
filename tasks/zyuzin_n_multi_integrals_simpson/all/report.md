# Вычисление многомерных интегралов с использованием многошаговой схемы (метод Симпсона) — ALL

- Student: Зюзин Никита Михайлович
- Technology: ALL (MPI + OpenMP)
- Variant: 11

## 1. Контекст

Гибридная версия объединяет два независимых уровня параллелизма:

- Межпроцессный уровень (MPI) — линейное пространство узлов сетки делится между процессами,
  каждый обрабатывает свой непрерывный диапазон;
- Внутрипроцессный уровень (OpenMP) — внутри каждого MPI-процесса запускается OMP-редукция
  по локальному поддиапазону.

Такой подход позволяет эффективно задействовать как несколько узлов кластера, так и несколько
ядер внутри каждого узла. Последовательная версия описана в `seq/report.md`, однопоточные
параллельные варианты — в `omp/report.md`, `stl/report.md`, `tbb/report.md`.

## 2. Постановка задачи

### Входные данные

Структура `SimpsonInput` из `common/include/common.hpp`:

- `lower_bounds` — вектор нижних пределов интегрирования;
- `upper_bounds` — вектор верхних пределов интегрирования;
- `n_steps` — вектор числа шагов (чётные, положительные);
- `func` — подынтегральная функция `std::function<double(const std::vector<double>&)>`.

### Выходные данные

Значение интеграла типа `double`.

### Ограничения

- Размеры всех трёх векторов должны совпадать и быть ненулевыми;
- `lower_bounds[i] <= upper_bounds[i]`;
- `n_steps[i] > 0` и `n_steps[i] % 2 == 0`;
- `func` должна быть непустой.

## 3. Базовый алгоритм

Математика идентична SEQ: составной многомерный метод Симпсона, `total_points = ∏(n_steps[dim]+1)`,
линейная нумерация, декодирование через смешанную систему счисления, тензорные веса Симпсона,
масштабирование `∏(h[dim]/3)`.

## 4. Межпроцессная схема

### Роли процессов

Все MPI-процессы равноправны: каждый самостоятельно вычисляет своё подмножество узлов и
сообщает локальную сумму. Нет выделенного «мастер»-процесса для раздачи данных — входные
данные одинаковы у всех рангов.

### Блочное разбиение по рангам

```cpp
int rank = 0, size = 1;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

const size_t usize = static_cast<size_t>(size);
const size_t urank = static_cast<size_t>(rank);
const size_t begin = (total_points * urank) / usize;
const size_t end   = (total_points * (urank + 1U)) / usize;
```

Ранг `r` из `P` обрабатывает диапазон `[N*r/P, N*(r+1)/P)`. Разница размеров диапазонов между
рангами не превышает 1.

### Используемые MPI-вызовы

| Вызов         | Место                      | Назначение                                             |
|---------------|----------------------------|--------------------------------------------------------|
| MPI_Comm_rank | RunImpl                    | Определение ранга текущего процесса                    |
| MPI_Comm_size | RunImpl                    | Определение общего числа процессов                     |
| MPI_Allreduce | RunImpl (после OMP-секции) | Суммирование локальных сумм всех рангов в глобальную   |

`MPI_Allreduce` выполняется симметрично на всех рангах — все получают итоговую сумму:

```cpp
double global_sum = 0.0;
MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```

Нет `MPI_Scatter`/`MPI_Gather` входных данных: каждый ранг сам вычисляет свои индексы узлов
из `total_points`, `rank` и `size`.

## 5. Внутрипроцессная схема

Функция `ComputeLocalSimpsonSum(begin, end)` параллелизует цикл по локальному диапазону ранга:

```cpp
#pragma omp parallel for default(none) shared(input, h, begin, end, num_dims) reduction(+ : local_sum)
for (size_t point_idx = begin; point_idx < end; ++point_idx) {
    auto temp = point_idx;
    std::vector<double> point(num_dims);
    double weight = 1.0;
    for (size_t dim = 0; dim < num_dims; ++dim) {
        const auto axis_points = static_cast<size_t>(input.n_steps[dim]) + 1U;
        const auto index = static_cast<int>(temp % axis_points);
        temp /= axis_points;
        point[dim] = input.lower_bounds[dim] + (static_cast<double>(index) * h[dim]);
        weight *= GetSimpsonWeight(index, input.n_steps[dim]);
    }
    local_sum += weight * input.func(point);
}
```

Параметры OpenMP:

- `default(none)` — запрет неявного захвата;
- `shared(input, h, begin, end, num_dims)` — только чтение;
- `reduction(+ : local_sum)` — локальная редукция без явной синхронизации.

### Финальное масштабирование

После `MPI_Allreduce` каждый ранг самостоятельно вычисляет масштабирующий фактор:

```cpp
double factor = ∏(h[dim] / 3.0);
result_ = global_sum * factor;
```

## 6. Детали реализации

Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`

### Методы класса `ZyuzinNSimpsonALL`

`ValidationImpl()` — проверка размеров, пределов, чётности шагов, непустоты `func`.

`PreProcessingImpl()` — сброс `result_` в 0.0.

`RunImpl()` — вычисляет `total_points`, определяет диапазон ранга, вызывает
`ComputeLocalSimpsonSum`, выполняет `MPI_Allreduce`, масштабирует и сохраняет результат.

`PostProcessingImpl()` — запись `result_` в `GetOutput()`.

`ComputeLocalSimpsonSum(size_t begin, size_t end)` — OMP-параллельный цикл по локальному
диапазону с редукцией.

`GetSimpsonWeight(int index, int n)` — статическая вспомогательная функция.

## 7. Проверка корректности

Те же 12 функциональных тестов с допуском `1e-3`. При запуске под `mpirun` тесты выполняются
всеми рангами, каждый обрабатывает свой поддиапазон, итоговый результат собирается через
`MPI_Allreduce`. Все 12 тестов пройдены при запуске с `mpirun -np 2`.

| ID | Название       | Функция        | Область              | Ожидание  |
|----|----------------|----------------|----------------------|-----------|
| 0  | 1d_linear      | x              | [0,1]                | 0.5       |
| 1  | 1d_quadratic   | x^2            | [0,1]                | 1/3       |
| 2  | 2d_sum         | x+y            | [0,1]^2              | 1.0       |
| 3  | 2d_product     | x·y            | [0,1]^2              | 0.25      |
| 4  | 2d_sum_squares | x^2+y^2        | [0,1]^2              | 2/3       |
| 5  | 2d_constant    | 1              | [0,2]×[0,3]          | 6.0       |
| 6  | 3d_constant    | 1              | [0,1]^3              | 1.0       |
| 7  | 3d_sum         | x+y+z          | [0,1]^3              | 1.5       |
| 8  | 1d_sin         | sin(x)         | [0,π]                | 2.0       |
| 9  | 2d_sin_cos     | sin(x)·cos(y)  | [0,π/2]^2            | 1.0       |
| 10 | 2d_exp         | exp(x+y)       | [0,1]^2              | (e−1)^2   |
| 11 | 3d_product     | x·y·z          | [0,1]^3              | 1/8       |

Performance-тест: `sin(x)·cos(y)·exp(z)` на `[0,π]×[0,π/2]×[0,1]`,
`n_steps = {180, 180, 120}`, ~4M узлов, ожидание `2·(e−1) ≈ 3.4366`.

## 8. Среда и результаты

Hardware:

- Процессор: AMD Ryzen 5 5600X
- Ядра/потоки: 6 ядер / 12 потоков
- ОЗУ: 16 GB
- ОС: Ubuntu 25.10
- Архитектура: x64

Toolchain:

- Компилятор: GCC 15.2.0
- IDE: Visual Studio Code 2026
- Тип сборки: Release
- Система сборки: CMake
- Версия MPI: Open MPI 5.0.8

Переменные окружения: `OMP_NUM_THREADS=N`, `mpirun -np P`

Результаты (Release-сборка, масштабирование `ALL` до 12 workers):

| Mode       | Процессы | OMP-потоков | Время, с | Ускорение | Эффективность |
|------------|----------|-------------|----------|-----------|---------------|
| seq (ref)  | 1        | 1           | 1.0635   | 1.00      | 100%          |
| all (task) | 2        | 1           | 0.3328   | 3.20      | 160%          |
| all (pipe) | 2        | 1           | 0.3243   | 3.30      | 165%          |
| all (task) | 2        | 2           | 0.2710   | 3.92      | 98%           |
| all (pipe) | 2        | 2           | 0.2711   | 3.94      | 99%           |
| all (task) | 2        | 3           | 0.2748   | 3.87      | 64%           |
| all (pipe) | 2        | 3           | 0.2712   | 3.94      | 66%           |
| all (task) | 2        | 4           | 0.2700   | 3.94      | 49%           |
| all (pipe) | 2        | 4           | 0.2698   | 3.96      | 50%           |
| all (task) | 2        | 5           | 0.2763   | 3.85      | 39%           |
| all (pipe) | 2        | 5           | 0.2759   | 3.87      | 39%           |
| all (task) | 2        | 6           | 0.2700   | 3.94      | 33%           |
| all (pipe) | 2        | 6           | 0.2696   | 3.96      | 33%           |
| all (task) | 4        | 1           | 0.1725   | 6.17      | 154%          |
| all (pipe) | 4        | 1           | 0.1745   | 6.12      | 153%          |
| all (task) | 4        | 2           | 0.1376   | 7.73      | 97%           |
| all (pipe) | 4        | 2           | 0.1381   | 7.74      | 97%           |
| all (task) | 4        | 3           | 0.1384   | 7.69      | 64%           |
| all (pipe) | 4        | 3           | 0.1381   | 7.74      | 64%           |

## 9. Репродуцируемость

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

## 10. Выводы

Гибридная ALL-версия корректна: все 12 функциональных тестов пройдены при `mpirun -np 2`.
Лучшая конфигурация в диапазоне до 12 workers — `4×2` (4 MPI-процесса,
2 OMP-потока), почти идентичный результат показывает `4×3`.

Сильные стороны:

- Два уровня параллелизма: MPI делит узлы между процессами, OMP — между потоками внутри процесса;
- Симметричная сборка через `MPI_Allreduce` — нет узкого места в виде ранга-мастера;
- Нет передачи входных данных между рангами: каждый ранг вычисляет свой диапазон независимо.

Слабые стороны:

- Каждый ранг хранит полную копию входных данных (нет `MPI_Scatter`);
- Per-iteration выделение `std::vector<double> point` внутри OMP-цикла (как в OMP-версии);
- Для малого числа узлов overhead MPI-инициализации может превышать выигрыш.

### Когда гибридная схема оправдана

- Большое число узлов (большие `n_steps`): линейное масштабирование по числу процессов;
- Многоузловые вычисления: каждый узел кластера запускает один MPI-процесс с OMP-потоками
  на все доступные ядра.

### Не оправдана, когда

- Задача малой размерности или с небольшим числом шагов: overhead MPI превышает выигрыш;
- Доступна только одна машина с несколькими ядрами: OMP или TBB дают лучшее соотношение
  производительности к сложности.
