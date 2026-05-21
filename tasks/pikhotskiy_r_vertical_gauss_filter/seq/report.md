# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — SEQ

- Student: Пихотский Роман Владимирович, group 3823Б1ФИ1
- Technology: SEQ
- Variant: 25

## 1. Introduction
В этой версии реализован базовый последовательный алгоритм фильтрации изображения ядром Гаусса 3x3.  
SEQ-реализация используется как эталон корректности для параллельных версий.

## 2. Problem Statement
- Вход: изображение в градациях серого, заданное как `width`, `height` и одномерный массив `data`.
- Выход: изображение того же размера после сглаживания Гауссом.
- Ограничения:
  - `width > 0`, `height > 0`;
  - `data.size() == width * height`.
- Границы обрабатываются через `clamp` к ближайшему валидному индексу.

## 3. Baseline Algorithm (Sequential)
Используется разложение 3x3 ядра Гаусса на два одномерных прохода:
1. Вертикальный проход по каждой колонке: коэффициенты `[1, 2, 1]`.
2. Горизонтальный проход по каждой строке: коэффициенты `[1, 2, 1]`.
3. Нормализация результата: `(sum + 15) / 16`.

Вертикальное разбиение в SEQ организовано полосами по оси `x`, но вычисления идут в одном потоке.

## 4. Parallelization Scheme
Параллелизма нет.  
Логика разбиения на вертикальные полосы сохранена, чтобы структура кода совпадала с параллельными реализациями.

## 5. Implementation Details
- Основной класс: `PikhotskiyRVerticalGaussFilterSEQ`.
- Файлы:
  - `seq/include/ops_seq.hpp`
  - `seq/src/ops_seq.cpp`
- Этапы пайплайна:
  - `ValidationImpl` — проверка размеров и размера массива;
  - `PreProcessingImpl` — подготовка буферов;
  - `RunImpl` — 2 прохода (vertical/horizontal);
  - `PostProcessingImpl` — сбор результата в `GetOutput()`.
- Для индексации используется линейный индекс `row * width + col`.

## 6. Experimental Setup
- CPU/OS/Compiler зависят от стенда запуска (локально или CI).
- Сборка в `Release`.
- Для проверок используются:
  - функциональные тесты: `tests/functional/main.cpp`;
  - performance-тесты: `tests/performance/main.cpp`.

Пример команд:
```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target ppc_func_tests ppc_perf_tests
```

## 7. Results and Discussion
### 7.1 Correctness
Корректность проверяется сравнением с заранее заданными ожидаемыми матрицами в функциональных тестах, включая:
- малые размеры (`1x1`, `2x2`, `3x3`);
- граничные случаи (`zero width/height`, `data size mismatch`).

### 7.2 Performance
SEQ служит baseline для сравнений с OMP/TBB/STL/ALL.  
Метрики скорости считаются в `tests/performance/main.cpp` через инфраструктуру курса.

## 8. Conclusions
SEQ-версия полностью решает задачу и формирует корректный эталон для всех параллельных реализаций.  
Двухпроходная схема уменьшает вычислительную сложность относительно прямой 3x3 свертки.

## 9. References
1. OpenCV Gaussian filtering concept: <https://docs.opencv.org/>
2. Parallel Programming Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
