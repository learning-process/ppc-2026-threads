# Вычисление многомерных интегралов методом прямоугольников

- Студент: Дергунов Сергей Антонович, группа 3823Б1ПР4
- Задача: `dergynov_s_integrals_multistep_rectangle`
- Вариант: 9
- Реализации: SEQ, OMP, TBB, STL, ALL

## 1. Введение

Реализован алгоритм вычисления многомерных интегралов методом средних
прямоугольников. Поддержаны пять технологий параллельного программирования:
`seq`, `omp`, `tbb`, `stl`, `all`.

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

Общий алгоритм для всех реализаций (метод средних прямоугольников):

1. Вычисляется шаг сетки по каждому измерению: `h_i = (right_i - left_i) / n`
2. Объём элементарной ячейки: `cell_volume = произведение h_i`
3. Перебираются все узлы сетки:
   - Для каждого узла вычисляются координаты средней точки ячейки
   - Вычисляется значение функции в этой точке
   - Сумма накапливается
4. Результат: `integral = sum * cell_volume`

Асимптотика: O(n^d) по времени, O(d) по памяти (где d — размерность).

## 4. Схема распараллеливания

| Технология | Схема распараллеливания                                                          |
|------------|----------------------------------------------------------------------------------|
| SEQ        | Без параллелизма, baseline                                                       |
| OMP        | `#pragma omp parallel for schedule(static)` по точкам сетки                      |
| TBB        | `tbb::parallel_reduce` с `tbb::blocked_range`                                    |
| STL        | `std::thread`, равномерное разбиение точек между потоками                        |
| ALL        | MPI + OpenMP: распределение точек между процессами, внутри процесса — по потокам |

## 5. Детали реализации

- **Общие типы**: `common/include/common.hpp`
- **Папки реализаций**: `seq`, `omp`, `tbb`, `stl`, `all`
- **Единый конвейер**: `ValidationImpl`, `PreProcessingImpl`, `RunImpl`, `PostProcessingImpl`

## 6. Экспериментальная среда

### Окружение

- **CPU**: AMD Ryzen 3 3200G
- **RAM**: 16 GB
- **OS**: Linux (Docker container)
- **Compiler**: Clang 21.1.8
- **Тип сборки**: Release

### Сборка

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build --parallel
```

### Запуск функциональных тестов

```bash
mpirun --allow-run-as-root -n 4 ./build/bin/ppc_func_tests --gtest_filter=*Dergynov*
```

### Запуск тестов производительности

```bash
export PPC_NUM_THREADS=1
./build/bin/ppc_perf_tests --gtest_filter="*seq*"

export PPC_NUM_THREADS=2
./build/bin/ppc_perf_tests --gtest_filter="*omp*"
./build/bin/ppc_perf_tests --gtest_filter="*tbb*"
./build/bin/ppc_perf_tests --gtest_filter="*stl*"
./build/bin/ppc_perf_tests --gtest_filter="*all*"

export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*omp*"
./build/bin/ppc_perf_tests --gtest_filter="*tbb*"
./build/bin/ppc_perf_tests --gtest_filter="*stl*"
./build/bin/ppc_perf_tests --gtest_filter="*all*"
```

## 7. Результаты

### 7.1 Корректность

Проверка корректности выполняется функциональными тестами.
Параллельные реализации сравниваются с baseline SEQ.
Все 95 тестов пройдены успешно.

### 7.2 Производительность

Определения метрик:

- `workers` — число исполнительных единиц (потоков)
- `time` — wall-clock время, секунды (режим task_run)
- `speedup = T_seq / T_mode`
- `efficiency = speedup / workers * 100%`

Сводная таблица результатов:

| mode | workers | time, s | speedup | efficiency, % |
|------|--------:|--------:|--------:|--------------:|
| seq  |       1 | 0.02584 | 1.00    | N/A           |
| omp  |       2 | 0.01805 | 1.43    | 71.5%         |
| omp  |       4 | 0.02037 | 1.27    | 31.7%         |
| tbb  |       2 | 0.02810 | 0.92    | 46.0%         |
| tbb  |       4 | 0.02127 | 1.21    | 30.3%         |
| stl  |       2 | 0.02671 | 0.97    | 48.5%         |
| stl  |       4 | 0.01530 | 1.69    | 42.2%         |
| all  |       2 | 0.03955 | 0.65    | 32.5%         |
| all  |       8 | 0.03676 | 0.70    | 8.8%          |

### 7.3 Интерпретация результатов

- **OMP**: Наибольшее ускорение при 2 потоках (1.43x). При 4 потоках ускорение падает из-за накладных расходов на синхронизацию.
- **TBB**: Ускорение нестабильно: 0.92x (2 потока), 1.21x (4 потока) из-за накладных расходов runtime.
- **STL**: Лучшее ускорение — 1.69x (4 потока) благодаря эффективному ручному управлению потоками.
- **ALL**: Гибридная MPI+OpenMP версия показывает ускорение 0.70x при 8 workers. Накладные расходы на
коммуникацию ограничивают масштабируемость.

## 8. Выводы

Задача реализована и документирована в едином стиле для всех пяти технологий.
SEQ используется как baseline. Среди параллельных версий наилучшее ускорение
показывает STL (1.69x при 4 потоках). OMP и TBB дают умеренное ускорение
(до 1.43x). ALL версия демонстрирует гибридную модель MPI+OpenMP.

## 9. Источники

1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. OpenMP specification: <https://www.openmp.org/specifications/>
3. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
4. C++ reference (std::thread): <https://en.cppreference.com/w/cpp/thread/thread>
5. Course report requirements: `docs/common_information/report.rst`
