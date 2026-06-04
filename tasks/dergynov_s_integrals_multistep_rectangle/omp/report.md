# Вычисление многомерных интегралов методом прямоугольников — OMP

- Студент: Дергунов Сергей Антонович, группа 3823Б1ПР4
- Технология: OMP
- Вариант: 9

## 1. Введение

OMP-версия использует OpenMP для распараллеливания цикла по точкам сетки.
Математика вычислений полностью совпадает с SEQ.

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

Алгоритм совпадает с SEQ (метод средних прямоугольников).

## 4. Схема распараллеливания

Используется `#pragma omp parallel for schedule(static)` для распределения
итераций по точкам сетки между потоками.

- **shared**: функция, границы, шаги, размерность, общее количество точек
- **private**: локальная сумма, координаты точки
- **reduction**: не используется (локальные суммы собираются в массив)
- **schedule(static)**: фиксированное распределение, подходит для равномерной нагрузки
- После `parallel for` действует неявный барьер

## 5. Детали реализации

- Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`
- Класс: `DergynovSIntegralsMultistepRectangleOMP`
- Конвейер: `ValidationImpl`, `PreProcessingImpl`, `RunImpl`, `PostProcessingImpl`

## 6. Экспериментальная среда

Сборка:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON
cmake --build build --parallel
```

Запуск:

```bash
export PPC_NUM_THREADS=2
export OMP_NUM_THREADS=2
./build/bin/ppc_perf_tests --gtest_filter="*omp*"

export PPC_NUM_THREADS=4
export OMP_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*omp*"
```

## 7. Результаты

### 7.1 Корректность

OMP-версия проходит те же функциональные тесты, что и baseline SEQ.
Результаты совпадают с точностью до 1e-6.

### 7.2 Производительность

Определения метрик:

- `workers` — число потоков OpenMP
- `time` — wall-clock время, секунды
- `speedup = T_seq / T_mode`
- `efficiency = speedup / workers * 100%`

Результаты:

| mode | workers | time (task_run), s | time (pipeline), s | speedup | efficiency, % |
|------|--------:|-------------------:|-------------------:|--------:|--------------:|
| seq  |       1 |            0.02584 |            0.02783 |    1.00 |           N/A |
| omp  |       2 |            0.01805 |            0.01972 |    1.43 |         71.5% |
| omp  |       4 |            0.02037 |            0.02882 |    1.27 |         31.7% |

Комментарий: Наибольшее ускорение достигается при 2 потоках (1.43x).
При 4 потоках ускорение падает из-за накладных расходов на синхронизацию
и ограниченной параллелимости задачи.

## 8. Выводы

OMP-версия сохраняет корректность baseline и обеспечивает ускорение
до 1.43x при 2 потоках.

## 9. Источники

1. OpenMP specification: <https://www.openmp.org/specifications/>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
