# Линейная фильтрация изображений (вертикальное разбиение). Ядро Гаусса 3x3

- Student: Пихотский Роман Владимирович, group 3823Б1ФИ1
- Task: `pikhotskiy_r_vertical_gauss_filter`
- Variant: 25
- Implementations: SEQ, OMP, TBB, STL, ALL

## 1. Introduction
Цель работы — реализовать фильтрацию изображения ядром Гаусса 3x3 с вертикальным разбиением и подготовить несколько backend-ов в инфраструктуре курса.

## 2. Problem Statement
- Входные данные:
  - `width`, `height`;
  - `data` (`std::vector<std::uint8_t>`) длиной `width * height`.
- Выход:
  - изображение того же размера после применения фильтра Гаусса.
- Основные требования:
  - строгая валидация входа;
  - корректная обработка границ изображения;
  - совместимость с `BaseTask` и тестовым пайплайном курса.

## 3. Baseline Algorithm (Sequential)
Для вычислений используется сепарабельное представление ядра Гаусса:

```text
[1 2 1]
[2 4 2] / 16
[1 2 1]
```

Вместо прямой 3x3 свертки выполняются два прохода:
1. Вертикальный `[1, 2, 1]` в промежуточный `int`-буфер.
2. Горизонтальный `[1, 2, 1]` в итоговый `uint8_t`-буфер.

Границы: `clamp` индексов в диапазон `[0, width-1]` / `[0, height-1]`.

## 4. Parallelization Scheme
- **SEQ**: однопоточная обработка полос.
- **OMP**: `parallel for` по вертикальным полосам, `schedule(static)`.
- **TBB**: `oneapi::tbb::parallel_for` + `blocked_range` по полосам.
- **STL**: `std::thread`, явное распределение полос между worker-потоками.
- **ALL**: backend `kALL`, использующий ту же схему полос и потоковый запуск проходов.

Во всех backend-ах сохранена одинаковая математическая модель, что упрощает проверку корректности.

## 5. Implementation Details
- Папка задачи: `tasks/pikhotskiy_r_vertical_gauss_filter`.
- Общий тип данных: `common/include/common.hpp` (`Matrix`).
- Реализованные backend-классы:
  - `PikhotskiyRVerticalGaussFilterSEQ`
  - `PikhotskiyRVerticalGaussFilterOMP`
  - `PikhotskiyRVerticalGaussFilterTBB`
  - `PikhotskiyRVerticalGaussFilterSTL`
  - `PikhotskiyRVerticalGaussFilterALL`
- Все версии следуют единому pipeline:
  - `ValidationImpl`
  - `PreProcessingImpl`
  - `RunImpl`
  - `PostProcessingImpl`

## 6. Experimental Setup
- Сборка: `Release`.
- Проверки:
  - `tests/functional/main.cpp`
  - `tests/performance/main.cpp`
- Управление конфигурацией:
  - `PPC_NUM_THREADS` для потоковых backend-ов.

Пример сборки:
```bash
cmake -S . -B build -DUSE_FUNC_TESTS=ON -DUSE_PERF_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target ppc_func_tests ppc_perf_tests
```

## 7. Results and Discussion
### 7.1 Correctness
Корректность подтверждается:
- параметризованными функциональными тестами на малых матрицах;
- негативными тестами на невалидный ввод (`zero sizes`, `data size mismatch`);
- сравнением выходов разных backend-ов с эталоном SEQ.

### 7.2 Performance
Performance-тесты включают все backend-ы (`SEQ/OMP/TBB/STL/ALL`) и позволяют сравнить время выполнения в одинаковых условиях.
Фактические численные значения зависят от аппаратуры и параметров запуска.

## 8. Conclusions
Задача реализована во всех требуемых технологиях с единым алгоритмическим ядром и едиными проверками корректности.  
Это обеспечивает сопоставимость результатов и упрощает сопровождение решения.

## 9. References
1. Parallel Programming Course repository: <https://github.com/learning-process/ppc-2026-threads>
2. OpenMP specification: <https://www.openmp.org/specifications/>
3. oneTBB documentation: <https://uxlfoundation.github.io/oneTBB/>
4. C++ reference (`std::thread`): <https://en.cppreference.com/w/cpp/thread/thread>
5. Report requirements: `docs/common_information/report.rst`
