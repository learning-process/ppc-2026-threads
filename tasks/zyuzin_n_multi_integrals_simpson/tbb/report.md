# Вычисление многомерных интегралов с использованием многошаговой схемы (метод Симпсона) — TBB

- Student: Зюзин Никита Михайлович
- Technology: TBB (oneTBB)
- Variant: 11

## 1. Контекст

TBB-версия использует `oneapi::tbb::parallel_reduce` для параллельной редукции по пространству
узлов интегральной сетки. По сравнению с OMP-версией предоставляет более тонкий контроль над
гранулярностью задач через `grain_size`. Последовательная версия описана в `seq/report.md`,
альтернативы — в `omp/report.md`, `stl/report.md`, `all/report.md`.

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

## 4. Схема параллелизма

### Контроль числа потоков

```cpp
const int max_threads = static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
oneapi::tbb::global_control thread_limiter(
    oneapi::tbb::global_control::max_allowed_parallelism, max_threads);
```

`global_control` ограничивает максимальный параллелизм TBB числом аппаратных потоков,
предотвращая избыточное создание задач.

### Параллельная редукция с grain_size

```cpp
constexpr size_t kGrainSize = 4096;
const double sum = oneapi::tbb::parallel_reduce(
    oneapi::tbb::blocked_range<size_t>(0, total_points, kGrainSize),
    0.0,
    [&](const oneapi::tbb::blocked_range<size_t> &range, double local_sum) {
        std::vector<double> point(num_dims);   // один раз на grain
        for (size_t point_idx = range.begin(); point_idx < range.end(); ++point_idx) {
            auto temp = point_idx;
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
        return local_sum;
    },
    [](double left, double right) { return left + right; }
);
```

### Ключевые параметры

- `kGrainSize = 4096` — каждая TBB-задача обрабатывает минимум 4096 узлов. Для задачи с
  ~4M узлов это даёт не более ~977 задач, что снижает накладные расходы TBB-планировщика до
  приемлемого уровня;
- `std::vector<double> point(num_dims)` выделяется **один раз на grain**, а не на итерацию
  (оптимизация относительно SEQ/OMP);
- Комбинирующая лямбда `[](double left, double right) { return left + right; }` — чистое
  суммирование без дополнительных аллокаций.

### Балансировка нагрузки

TBB использует work-stealing планировщик: если один поток завершил свои grain-задачи раньше
других, он «крадёт» задачи из очередей других потоков. Для данной задачи (равномерная нагрузка
на узел) статическое распределение и work-stealing дают сопоставимые результаты.

## 5. Детали реализации

Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`

### Методы класса `ZyuzinNSimpsonTBB`

`ValidationImpl()` — проверка размеров, пределов, чётности шагов, непустоты `func`.

`PreProcessingImpl()` — сброс `result_` в 0.0.

`RunImpl()` — вызов `ComputeSimpsonMultiDim()`.

`PostProcessingImpl()` — запись `result_` в `GetOutput()`.

`ComputeSimpsonMultiDim()` — устанавливает `global_control`, запускает `parallel_reduce`
по диапазону `[0, total_points)`.

`GetSimpsonWeight(int index, int n)` — статическая вспомогательная функция.

### Зависимости

```cpp
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_reduce.h"
```

## 6. Проверка корректности

12 функциональных тестов с допуском `1e-3`. Корректность `parallel_reduce` обеспечивается
ассоциативностью и коммутативностью суммирования. Порядок суммирования grain-результатов может
отличаться от SEQ, поэтому допуск `1e-3` (а не строгое равенство) достаточен для всех тестов.

| ID | Название       | Функция        | Область              | Ожидание  |
|----|----------------|----------------|----------------------|-----------|
| 0  | 1d_linear      | x              | [0,1]                | 0.5       |
| 1  | 1d_quadratic   | x²             | [0,1]                | 1/3       |
| 2  | 2d_sum         | x+y            | [0,1]²               | 1.0       |
| 3  | 2d_product     | x·y            | [0,1]²               | 0.25      |
| 4  | 2d_sum_squares | x²+y²          | [0,1]²               | 2/3       |
| 5  | 2d_constant    | 1              | [0,2]×[0,3]          | 6.0       |
| 6  | 3d_constant    | 1              | [0,1]³               | 1.0       |
| 7  | 3d_sum         | x+y+z          | [0,1]³               | 1.5       |
| 8  | 1d_sin         | sin(x)         | [0,π]                | 2.0       |
| 9  | 2d_sin_cos     | sin(x)·cos(y)  | [0,π/2]²             | 1.0       |
| 10 | 2d_exp         | exp(x+y)       | [0,1]²               | (e−1)²   |
| 11 | 3d_product     | x·y·z          | [0,1]³               | 1/8       |

Performance-тест: `sin(x)·cos(y)·exp(z)` на `[0,π]×[0,π/2]×[0,1]`,
`n_steps = {180, 180, 120}`, ~4M узлов, ожидание `2·(e−1) ≈ 3.4366`.

## 7. Среда и результаты

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

Переменная окружения: `PPC_NUM_THREADS=N`

| Mode       | Потоки | Время, с | Ускорение | Эффективность |
|------------|--------|----------|-----------|---------------|
| seq (task) | 1      | 1.1135   | 1.00      | 100%          |
| tbb (task) | 2      | 0.1972   | 5.65      | 283%          |
| tbb (pipe) | 2      | 0.1987   | 5.40      | 270%          |
| tbb (task) | 4      | 0.1028   | 10.83     | 271%          |
| tbb (pipe) | 4      | 0.1045   | 10.27     | 257%          |
| tbb (task) | 6      | 0.0721   | 15.45     | 258%          |
| tbb (pipe) | 6      | 0.0738   | 14.54     | 242%          |
| tbb (task) | 8      | 0.0728   | 15.30     | 191%          |
| tbb (pipe) | 8      | 0.0722   | 14.87     | 186%          |
| tbb (task) | 10     | 0.0698   | 15.94     | 159%          |
| tbb (pipe) | 10     | 0.0686   | 15.65     | 156%          |
| tbb (task) | 12     | 0.0677   | 16.46     | 137%          |
| tbb (pipe) | 12     | 0.0667   | 16.10     | 134%          |

## 8. Репродуцируемость

Команды запуска функциональных тестов

```cpp
// SEQ
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*seq*"

// TBB (2 потока)
$env:PPC_NUM_THREADS=2
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*tbb*"
```

Команды запуска тестов производительности

```cpp
// SEQ (baseline)
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*seq*"

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
```

## 9. Выводы

TBB-реализация использует `parallel_reduce` с `grain_size = 4096`, что обеспечивает
эффективное разбиение ~4M узлов на задачи. Переиспользование буфера `point` на grain (а не
на итерацию) снижает нагрузку на аллокатор.

Результаты подтверждают устойчивое масштабирование до 12 потоков: минимальное время в серии
достигается при 12 потоках, а участок 6–12 потоков остаётся стабильно близким к лучшему.

Сильные стороны:

- `grain_size = 4096` хорошо амортизирует накладные расходы TBB-планировщика на задачах ~4M узлов;
- `parallel_reduce` с work-stealing обеспечивает автоматическую балансировку нагрузки;
- Буфер `point` переиспользуется в рамках одного grain.

Слабые стороны:

- Выбор `grain_size` подобран эмпирически; для других размерностей или числа шагов оптимальное
  значение может отличаться;
- Несколько более сложная структура кода по сравнению с OMP-версией.
