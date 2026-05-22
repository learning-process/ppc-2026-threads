# Вычисление многомерных интегралов методом прямоугольников — ALL

- Студент: Дергунов Сергей Антонович, группа 3823Б1ПР4
- Технология: ALL
- Вариант: 9

## 1. Введение

ALL-версия объединяет все четыре реализации (SEQ, OMP, TBB, STL)
в одном backend-е. Для потоковой конфигурации используется схема
`ranks × threads = 1 × N` (только внутрипроцессный параллелизм).

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

ALL-версия последовательно вызывает все четыре реализации и проверяет
совпадение их результатов с заданной точностью. При несовпадении
возвращается ошибка. Итоговый результат берётся из TBB-реализации
как наиболее эффективной.

Для `threads`-задачи используется конфигурация `ranks × threads = 1 × N`:
- `ranks = 1`: межпроцессного обмена нет
- `threads = N`: обработка идёт потоками

## 5. Детали реализации

- Файлы: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`
- Класс: `DergynovSIntegralsMultistepRectangleALL`
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
./build/bin/ppc_perf_tests --gtest_filter="*all*"

export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*all*"
```

## 7. Результаты

### 7.1 Корректность

ALL-версия проходит те же функциональные тесты, что и baseline SEQ.
Все четыре реализации дают совпадающие результаты.

### 7.2 Производительность

Определения метрик:
- `workers` — число потоков
- `time` — wall-clock время, секунды
- `speedup = T_seq / T_mode`
- `efficiency = speedup / workers * 100%`

Результаты:
| mode | workers | time (task_run), s | time (pipeline), s | speedup | efficiency, % |
|------|--------:|-------------------:|-------------------:|--------:|--------------:|
| seq  |       1 |            0.02584 |            0.02783 |    1.00 |           N/A |
| all  |       2 |            0.13077 |            0.12882 |    0.20 |         10.0% |
| all  |       4 |            0.11448 |            0.12468 |    0.23 |          5.8% |

Комментарий: ALL-версия показывает значительно более низкую
производительность из-за того, что последовательно вызывает все четыре
реализации. Это демонстрационный режим, не предназначенный для
высокой производительности.

## 8. Выводы

ALL-версия поддерживает единый запуск в инфраструктуре курса,
демонстрирует все четыре технологии и проверяет согласованность их
результатов.

## 9. Источники

1. Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. OpenMP specification: <https://www.openmp.org/specifications/>
3. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
4. C++ reference (std::thread): <https://en.cppreference.com/w/cpp/thread/thread>
5. Course report requirements: `docs/common_information/report.rst`