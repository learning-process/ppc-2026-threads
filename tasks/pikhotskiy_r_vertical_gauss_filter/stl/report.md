# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — STL

- Student: Пихотский Роман Владимирович, group 3823Б1ФИ1
- Technology: STL (`std::thread`)
- Variant: 25

## 1. Introduction
STL-версия использует явное управление потоками через `std::thread` для параллельной обработки изображения по вертикальным полосам.

## 2. Problem Statement
- Вход: `width`, `height`, `data` (`uint8_t`).
- Выход: сглаженное изображение той же размерности.
- Проверки:
  - размеры положительные;
  - размер массива совпадает с `width * height`.

## 3. Baseline Algorithm (Sequential)
Алгоритм эквивалентен SEQ:
1. Вертикальная свертка `[1, 2, 1]`.
2. Горизонтальная свертка `[1, 2, 1]`.
3. Нормализация `(sum + 15) / 16`.
4. Границы через `clamp`.

## 4. Parallelization Scheme
Схема распараллеливания:
- изображение разбивается на `stripe_count` полос;
- полосы распределяются между `actual_threads`;
- для каждого прохода запускается группа `std::thread`, затем `join()`.

В реализации используется вспомогательная функция `RunPassInParallel(...)`, чтобы избежать дублирования логики распределения полос.

## 5. Implementation Details
- Класс: `PikhotskiyRVerticalGaussFilterSTL`.
- Файлы:
  - `stl/include/ops_stl.hpp`
  - `stl/src/ops_stl.cpp`
- Особенности:
  - `actual_threads = min(requested_threads, stripe_count)`;
  - запись в буферы потокобезопасна, так как полосы не пересекаются.

## 6. Experimental Setup
- Сборка: `Release`.
- Управление числом потоков: `PPC_NUM_THREADS`.
- Проверки:
  - функциональные тесты на корректность;
  - performance-тесты для измерения времени.

Пример:
```bash
export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

## 7. Results and Discussion
### 7.1 Correctness
Проверка проходит теми же тест-кейсами, что и для SEQ/OMP/TBB, включая негативные случаи валидации.

### 7.2 Performance
STL-версия масштабируется за счет параллельной обработки полос, но чувствительна к:
- накладным расходам на создание потоков;
- размеру изображения и доле синхронизации между проходами.

## 8. Conclusions
STL-реализация обеспечивает переносимый и прозрачный контроль над потоками и сохраняет точность SEQ.  
Вариант хорошо подходит как базовый многопоточный backend без внешнего runtime.

## 9. References
1. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
2. Parallel Programming Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
