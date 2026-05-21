# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3 — TBB

- Student: Пихотский Роман Владимирович, group 3823Б1ФИ1
- Technology: TBB
- Variant: 25

## 1. Introduction
TBB-версия реализует ту же фильтрацию Гаусса 3x3, что и SEQ/OMP, но с использованием oneTBB (`parallel_for`).  
Цель — получить масштабирование на многопоточном CPU при сохранении корректности.

## 2. Problem Statement
- Вход: матрица яркостей `width x height` в виде одномерного массива `uint8_t`.
- Выход: сглаженное изображение того же размера.
- Требования к входу:
  - `width > 0`, `height > 0`;
  - `data.size() == width * height`.

## 3. Baseline Algorithm (Sequential)
Применяется сепарабельная свертка:
1. Вертикальный проход `[1, 2, 1]` в промежуточный буфер `int`.
2. Горизонтальный проход `[1, 2, 1]`.
3. Нормировка `(sum + 15) / 16`.
4. Границы: `clamp` индексов.

## 4. Parallelization Scheme
Параллелизация выполнена по вертикальным полосам:
- `stripe_width = max(1, width / 8)`;
- `stripe_count = ceil(width / stripe_width)`;
- `oneapi::tbb::parallel_for(blocked_range<int>(0, stripe_count), ...)`.

Алгоритм выполняет два независимых параллельных этапа:
- vertical-pass по полосам;
- horizontal-pass по полосам.

## 5. Implementation Details
- Класс: `PikhotskiyRVerticalGaussFilterTBB`.
- Файлы:
  - `tbb/include/ops_tbb.hpp`
  - `tbb/src/ops_tbb.cpp`
- Важные детали:
  - записи в `vertical_buffer_` и `result_buffer_` непересекающиеся между полосами;
  - для CI-совместимости соблюдены прямые include и стиль clang-tidy.

## 6. Experimental Setup
- Сборка: `Release`.
- Количество потоков TBB управляется runtime и настройками окружения.
- Для сравнения используются:
  - functional-тесты;
  - performance-тесты из `tests/performance/main.cpp`.

Пример:
```bash
export PPC_NUM_THREADS=4
./build/bin/ppc_perf_tests --gtest_filter="*pikhotskiy_r_vertical_gauss_filter*"
```

## 7. Results and Discussion
### 7.1 Correctness
Результат TBB сверяется с ожидаемыми матрицами в функциональных тестах и должен совпадать с SEQ.

### 7.2 Performance
При росте размера изображения TBB обычно выигрывает за счет распределения полос между worker-потоками.  
Реальный speedup зависит от размера изображения и memory bandwidth.

## 8. Conclusions
TBB-реализация корректно повторяет математическую модель SEQ и дает переносимый способ многопоточной обработки.  
Подход с `blocked_range` удобно масштабируется и хорошо сочетается с вертикальным разбиением.

## 9. References
1. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
2. Parallel Programming Course repository: <https://github.com/learning-process/ppc-2026-threads>
3. Course report requirements: `docs/common_information/report.rst`
