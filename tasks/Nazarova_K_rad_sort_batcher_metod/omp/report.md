# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) — OMP

- Student: Назарова Ксения Олеговна
- Technology: OMP
- Variant: 9

## 1. Контекст

OMP переносит основной цикл обхода сетки из последовательной реализации в параллельную область. Независимость ячеек
позволяет суммировать вклады параллельно, а общий результат получать через редукцию.

## 2. Постановка задачи

Постановка совпадает с `seq/report.md`: требуется вычислить приближённое значение многомерного интеграла методом
средних прямоугольников. Последовательная версия используется как эталон корректности и как знаменатель для расчёта
ускорения.

## 3. Базовый алгоритм

В отличие от SEQ, многомерный индекс не увеличивается счётчиком. В OMP заранее вычисляется `total_cells_`, после чего
каждая итерация линейного диапазона `[0, total_cells_)` преобразуется в набор координат по осям.

## 4. Схема распараллеливания

Параллелится внешний линейный диапазон ячеек. Каждая итерация независима: она строит собственную точку в центре ячейки
и добавляет значение функции в локальную копию суммы. Общая переменная `sum` объединяется директивой
`reduction(+ : sum)`.

Ключевой фрагмент из `omp/src/ops_omp.cpp`:

```cpp
#pragma omp parallel default(none) shared(input, step_sizes) firstprivate(dimension, total_cells) reduction(+ : sum)
{
  std::vector<double> point(dimension, 0.0);

#pragma omp for schedule(static)
  for (std::size_t linear_index = 0U; linear_index < total_cells; ++linear_index) {
    std::size_t current_index = linear_index;
    for (std::size_t axis = 0U; axis < dimension; ++axis) {
      const std::size_t coordinate_index = current_index % input.steps[axis];
      current_index /= input.steps[axis];
      point[axis] = input.lower_bounds[axis] + ((static_cast<double>(coordinate_index) + 0.5) * step_sizes[axis]);
    }
    sum += input.function(point);
  }
}
```

`default(none)` заставляет явно указать область видимости переменных. `input` и `step_sizes` доступны только для
чтения, поэтому являются `shared`. `dimension` и `total_cells` передаются как `firstprivate`. Вектор `point` создаётся
внутри параллельной области, значит у каждого потока есть собственный буфер координат. `schedule(static)` выбран из-за
равномерной стоимости итераций на тестовом наборе.

## 5. Детали реализации

Файлы реализации: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.

Валидация повторяет общие проверки и дополнительно защищает вычисление `total_cells` от переполнения `std::size_t`.
В `PreProcessingImpl` вычисляются шаги, объём ячейки и общее число ячеек. Синхронизация вручную не требуется:
суммирование защищено OpenMP-редукцией, а в конце `parallel`-области есть неявный барьер.

## 6. Проверка корректности

OMP-версия прошла 7 функциональных тестов из общего набора. Результаты сравнивались с аналитически ожидаемыми
значениями с заданным `eps`. Отдельных расхождений при `PPC_NUM_THREADS=4` не обнаружено.

## 7. Экспериментальная среда

Окружение: Linux 6.6.114.1-microsoft-standard-WSL2, Intel Core i5-1235U, 12 логических CPU, GCC 13.3.0. Основная
конфигурация: `PPC_NUM_THREADS=4`, `OMP_NUM_THREADS=4`.

Команды запуска:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Release

cmake --build build -j --parallel

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mpirun --oversubscribe -n 4 ./build/bin/ppc_perf_tests --gtest_filter='*nazarova_k_calc_integ_rectangles_omp*'
```

## 8. Результаты

| Mode     | Threads | Time, s      | Speedup | Efficiency |
| -------- | ------- | ------------ | ------- | ---------- |
| pipeline | 4       | 1.7958475872 | 0.41    | 0.10       |
| task_run | 4       | 1.7755917146 | 0.39    | 0.10       |

На выбранном тесте ускорение меньше единицы. Причина — лёгкое тело цикла: преобразование индекса и вычисление линейной
функции не компенсируют накладные расходы OpenMP-области и редукции.

## 9. Выводы

OpenMP-реализация корректно распараллеливает независимые ячейки и не содержит общей записи без редукции. Для более
тяжёлой функции или большего числа ячеек схема должна масштабироваться лучше, но в данном performance-тесте overhead
превышает выигрыш.
