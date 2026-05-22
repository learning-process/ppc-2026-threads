# Вычисление многомерных интегралов методом прямоугольников — TBB

- Студент: Дергунов Сергей Антонович, группа 3823Б1ПР4
- Технология: TBB
- Вариант: 9

## 1. Введение

TBB-версия использует task-based модель oneTBB для распараллеливания
вычислений через `tbb::parallel_reduce`.

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

Используется `tbb::parallel_reduce` с `tbb::blocked_range`.

- **blocked_range**: задаёт диапазон линейных индексов точек сетки
- **grainsize**: не задан явно, используется значение по умолчанию
- **partitioner**: стандартный `auto_partitioner`
- **reduction**: встроенная редукция суммирования

Первая фаза (`parallel_reduce`) полностью завершается до получения результата,
поэтому синхронизация не требуется.

## 5. Детали реализации

- Файлы: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`
- Класс: `DergynovSIntegralsMultistepRectangleTBB`
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
./build/bin/ppc_perf_tests --gtest_filter="*tbb*"

export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*tbb*"
```

## 7. Результаты

### 7.1 Корректность

TBB-версия проходит те же функциональные тесты, что и baseline SEQ.
Результаты совпадают с точностью до 1e-6.

### 7.2 Производительность

Определения метрик:
- `workers` — число рабочих потоков TBB runtime
- `time` — wall-clock время, секунды
- `speedup = T_seq / T_mode`
- `efficiency = speedup / workers * 100%`

Результаты:
| mode | workers | time (task_run), s | time (pipeline), s | speedup | efficiency, % |
|------|--------:|-------------------:|-------------------:|--------:|--------------:|
| seq  |       1 |            0.02584 |            0.02783 |    1.00 |           N/A |
| tbb  |       2 |            0.02810 |            0.03855 |    0.92 |         46.0% |
| tbb  |       4 |            0.02127 |            0.02067 |    1.21 |         30.3% |

Комментарий: TBB показывает нестабильное ускорение.
При 2 потоках ускорение ниже 1x, что может быть связано с накладными расходами
runtime. При 4 потоках наблюдается ускорение 1.21x.

## 8. Выводы

TBB-версия сохраняет корректность baseline.
Ускорение ограничено особенностями задачи и накладными расходами runtime.

## 9. Источники

1. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
2. Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`