# Вычисление многомерных интегралов с использованием многошаговой схемы (метод Симпсона) — STL

- Student: Зюзин Никита Михайлович
- Technology: STL (std::thread)
- Variant: 11

## 1. Контекст

STL-версия реализует параллелизм вручную через `std::thread` без внешних библиотек. Ключевое
отличие от OMP-версии — явное управление блочным разбиением и переиспользование буферов на поток.
Последовательная версия описана в `seq/report.md`, альтернативы — в `omp/report.md`,
`tbb/report.md`, `all/report.md`.

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

### Определение числа потоков

```cpp
const size_t requested_threads = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
const size_t num_threads = std::min(requested_threads, total_points);
```

Число потоков адаптируется к размеру задачи: если `total_points` меньше запрошенного числа
потоков, используются только `total_points` потоков, чтобы не порождать пустые потоки.

### Блочное разбиение

Каждый поток `i` обрабатывает непрерывный диапазон узлов:

```cpp
const size_t begin = (thread_id * total_points) / num_threads;
const size_t end   = ((thread_id + 1) * total_points) / num_threads;
```

Такое деление гарантирует, что размер чанков отличается не более чем на 1, и не создаёт
накладных расходов на динамическое планирование.

### Оптимизация буферов

В отличие от SEQ/OMP, `std::vector<double> point(num_dims)` выделяется **один раз на поток**
до начала внутреннего цикла. Индексы вычисляются в локальных переменных без хранения в отдельном
векторе:

```cpp
auto compute_chunk = [&](size_t thread_id) {
    const size_t begin = (thread_id * total_points) / num_threads;
    const size_t end   = ((thread_id + 1) * total_points) / num_threads;
    std::vector<double> point(num_dims);   // один раз на поток
    double local_sum = 0.0;

    for (size_t point_idx = begin; point_idx < end; ++point_idx) {
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
    partial_sums[thread_id] = local_sum;
};
```

### Запуск потоков

Поток 0 выполняется в вызывающем потоке без порождения нового `std::thread`:

```cpp
for (size_t thread_id = 1; thread_id < num_threads; ++thread_id) {
    workers.emplace_back(compute_chunk, thread_id);
}
compute_chunk(0);
for (auto &worker : workers) {
    worker.join();
}
```

### Агрегация

```cpp
const double sum = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
```

## 5. Детали реализации

Файлы: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`

### Методы класса `ZyuzinNSimpsonSTL`

`ValidationImpl()` — проверка размеров, пределов, чётности шагов, непустоты `func`.

`PreProcessingImpl()` — сброс `result_` в 0.0.

`RunImpl()` — вызов `ComputeSimpsonMultiDim()`.

`PostProcessingImpl()` — запись результата в `GetOutput()`.

`ComputeSimpsonMultiDim()` — создаёт `partial_sums` и `workers`, запускает лямбду
`compute_chunk` на каждом потоке, собирает результат через `std::accumulate`.

`GetSimpsonWeight(int index, int n)` — статическая вспомогательная функция.

## 6. Проверка корректности

12 функциональных тестов с допуском `1e-3`, идентичных SEQ и OMP. Отсутствие гонок
гарантируется тем, что каждый поток пишет исключительно в свой элемент `partial_sums[thread_id]`
и читает только разделяемые входные данные, которые на момент старта потоков не изменяются.

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
| stl (task) | 2      | 0.2033   | 5.48      | 274%          |
| stl (pipe) | 2      | 0.2297   | 4.67      | 234%          |
| stl (task) | 4      | 0.1145   | 9.73      | 243%          |
| stl (pipe) | 4      | 0.1062   | 10.11     | 253%          |
| stl (task) | 6      | 0.0747   | 14.91     | 248%          |
| stl (pipe) | 6      | 0.0781   | 13.75     | 229%          |
| stl (task) | 8      | 0.0882   | 12.63     | 158%          |
| stl (pipe) | 8      | 0.0856   | 12.54     | 157%          |
| stl (task) | 10     | 0.0807   | 13.79     | 138%          |
| stl (pipe) | 10     | 0.0784   | 13.69     | 137%          |
| stl (task) | 12     | 0.0752   | 14.82     | 124%          |
| stl (pipe) | 12     | 0.0753   | 14.25     | 119%          |

## 8. Репродуцируемость

Команды запуска функциональных тестов

```cpp
// SEQ
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*seq*"

// STL (2 потока)
$env:PPC_NUM_THREADS=2
./build/bin/ppc_func_tests --gtest_filter="*zyuzin_n*stl*"
```

Команды запуска тестов производительности

```cpp
// SEQ (baseline)
./build/bin/ppc_perf_tests --gtest_filter="*zyuzin_n*seq*"

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
```

## 9. Выводы

STL-реализация полностью корректна. Ключевая оптимизация — переиспользование буфера `point`
на поток: в отличие от SEQ/OMP, не выполняется `new`/`delete` на каждой из ~4M итераций.

Масштабирование не монотонно: лучшее время достигается при 6–12 потоках,
а при 8 потоках наблюдается локальная деградация. Это согласуется с чувствительностью
ручного статического разбиения к влиянию планировщика ОС и кэш-эффектам.

Сильные стороны:

- Нет зависимости от OpenMP или TBB — только стандартная библиотека C++;
- Переиспользование `point`-буфера снижает давление на аллокатор;
- Явный контроль числа потоков через `ppc::util::GetNumThreads()`.

Слабые стороны:

- Статическое блочное разбиение: если итерации неравнозначны по вычислительной стоимости,
  возникает дисбаланс нагрузки (для данной задачи нагрузка равномерна);
- Ручное управление потоками более многословно по сравнению с OMP/TBB.
