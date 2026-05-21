# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — OMP

- Student: Пихотский Роман Владимирович, group 3823Б1ФИ1
- Technology: OMP
- Variant: 25

## 1. Introduction
В OMP-версии реализована параллельная фильтрация Гауссом 3x3 с вертикальным разбиением изображения.  
Алгоритм полностью эквивалентен SEQ-версии по математике и формату результата.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`uint8_t`).
- Выход: отфильтрованное изображение той же размерности.
- Условия корректности:
  - размеры положительные;
  - `data.size() == width * height`.

## 3. Baseline Algorithm (Sequential)
Базовая схема совпадает с SEQ:
1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Обработка границ через `clamp`.

## 4. Parallelization Scheme
Распараллеливание выполняется по вертикальным полосам:
- изображение делится на `stripe_count` полос по столбцам;
- каждая полоса обрабатывается независимо;
- сначала параллельно запускается vertical-pass, затем horizontal-pass.

Используемые директивы:
- `#pragma omp parallel for default(none) schedule(static)`.

Это устраняет гонки, так как каждый поток пишет в непересекающуюся область буфера.

## 5. Implementation Details
- Класс: `PikhotskiyRVerticalGaussFilterOMP`.
- Файлы:
  - `omp/include/ops_omp.hpp`
  - `omp/src/ops_omp.cpp`
- `PreProcessingImpl`:
  - вычисляет ширину полосы из `ppc::util::GetNumThreads()`;
  - выделяет `vertical_buffer_` и `result_buffer_`.
- `RunImpl`:
  - проверка размеров буферов;
  - 2 параллельных прохода.

## 6. Experimental Setup
- Сборка в `Release`.
- Для ограничения числа потоков используются:
  - `PPC_NUM_THREADS`;
  - при необходимости `OMP_NUM_THREADS`.
- Проверки:
  - `tests/functional/main.cpp`;
  - `tests/performance/main.cpp`.

Пример запуска:
```bash
export PPC_NUM_THREADS=4
export OMP_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

## 7. Results and Discussion
### 7.1 Correctness
OMP-версия проверяется на тех же наборах данных, что и SEQ, и должна выдавать побайтно идентичный результат.

### 7.2 Performance
Ожидается ускорение на средних и больших изображениях благодаря независимой обработке полос.  
Предел масштабируемости определяется пропускной способностью памяти и накладными расходами OpenMP runtime.

## 8. Conclusions
OMP-реализация сохраняет точность SEQ и добавляет параллелизм без изменения интерфейса задачи.  
Вертикальное разбиение подходит для простой и безопасной декомпозиции работы между потоками.

## 9. References
1. OpenMP API Specification: <https://www.openmp.org/specifications/>
2. Parallel Programming Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
