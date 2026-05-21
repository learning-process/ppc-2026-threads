# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — ALL

- Student: Пихотский Роман Владимирович, group 3823Б1ФИ1
- Technology: ALL
- Variant: 25

## 1. Introduction
ALL-версия объединяет реализацию задачи в backend с типом `kALL`, сохраняя ту же вычислительную схему, что и в остальных вариантах.

## 2. Problem Statement
- Вход: изображение в формате (`width`, `height`, `data`).
- Выход: изображение того же размера после фильтра Гаусса 3x3.
- Валидность входа:
  - положительные размеры;
  - совпадение `data.size()` и `width * height`.

## 3. Baseline Algorithm (Sequential)
Выполняется сепарабельная свертка:
1. Вертикальный проход `[1, 2, 1]`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Округляющая нормировка `(sum + 15) / 16`.
4. Граничные условия через `clamp`.

## 4. Parallelization Scheme
В текущем `ALL`-решении используется потоковая декомпозиция `std::thread` по вертикальным полосам:
- расчет числа полос;
- равномерное распределение полос между workers;
- два последовательных этапа (vertical-pass, horizontal-pass), внутри каждого — параллельная обработка.

Для запуска проходов применяется вспомогательная функция `RunPassInParallel(...)`.

## 5. Implementation Details
- Класс: `PikhotskiyRVerticalGaussFilterALL`.
- Файлы:
  - `all/include/ops_all.hpp`
  - `all/src/ops_all.cpp`
- Архитектура полностью совместима с инфраструктурой курса:
  - `ValidationImpl`;
  - `PreProcessingImpl`;
  - `RunImpl`;
  - `PostProcessingImpl`.

## 6. Experimental Setup
- Сборка в `Release`.
- Управление потоками: `PPC_NUM_THREADS`.
- Проверки:
  - functional: проверка точности и валидации;
  - performance: сравнение с другими backend-ами.

Пример:
```bash
export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

## 7. Results and Discussion
### 7.1 Correctness
ALL-версия проверяется теми же функциональными тестами и должна выдавать результат, идентичный SEQ.

### 7.2 Performance
Производительность зависит от числа полос, размера изображения и затрат на создание/синхронизацию потоков.
Для небольших изображений overhead может быть сопоставим с полезной работой.

## 8. Conclusions
ALL-реализация корректно встроена в task pipeline и совместима с тестовой инфраструктурой.  
Поведение по корректности соответствует другим backend-ам задачи.

## 9. References
1. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
2. Parallel Programming Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
