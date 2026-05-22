# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) — SEQ

- Student: Назарова Ксения Олеговна
- Technology: SEQ
- Variant: 9

## 1. Контекст

Последовательная версия является эталоном для всех остальных backend-ов. Она задаёт общую постановку задачи, правила валидации входа, порядок обхода многомерной сетки и baseline-время для расчёта ускорения.

## 2. Постановка задачи

Нужно вычислить приближённое значение многомерного интеграла методом средних прямоугольников. Вход содержит функцию, нижние и верхние границы по каждой координате и число шагов по каждой оси. Выходом является одно значение `double`.

Вход считается корректным, если функция задана, размерность больше нуля, размеры `lower_bounds`, `upper_bounds` и `steps` совпадают, все границы конечны, число шагов положительно, а нижняя граница не превосходит верхнюю.

## 3. Базовый алгоритм

В `PreProcessingImpl` вычисляются шаги сетки и объём одной ячейки. В `RunImpl` создаётся вектор индексов многомерной сетки. На каждой итерации индексы преобразуются в координаты центра ячейки, значение функции добавляется к сумме, после чего индекс увеличивается как многомерный счётчик.

Временная сложность равна `O(D * N)`, где `D` — размерность, а `N = product(steps[i])` — общее число ячеек. Память равна `O(D)`, так как хранятся только текущая точка, текущие индексы и шаги по осям.

Ключевой фрагмент находится в `seq/src/ops_seq.cpp`:

```cpp
while (true) {
  for (std::size_t i = 0; i < dimension_; ++i) {
    point[i] = input.lower_bounds[i] + ((static_cast<double>(indices[i]) + 0.5) * step_sizes_[i]);
  }
  sum += input.function(point);

  std::size_t axis = 0U;
  while (axis < dimension_) {
    ++indices[axis];
    if (indices[axis] < input.steps[axis]) {
      break;
    }
    indices[axis] = 0U;
    ++axis;
  }
  if (axis == dimension_) {
    break;
  }
}
```

Этот цикл посещает каждую ячейку ровно один раз. Итоговая сумма умножается на `cell_volume_`, поэтому каждая ячейка вносит вклад `f(center) * volume`.

## 4. Детали реализации

Файлы реализации: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`.

`ValidationImpl` делегирует проверку функции `HasValidInput`. `PreProcessingImpl` заполняет `step_sizes_`, вычисляет `cell_volume_`, обнуляет результат и выходное значение. `RunImpl` выполняет полный последовательный обход сетки. `PostProcessingImpl` не требует дополнительных действий, потому что результат записывается в `GetOutput()` сразу после вычисления.

## 5. Проверка корректности

Корректность проверялась общими функциональными тестами. Они покрывают константную функцию, линейную функцию, произведение координат, сумму координат, квадрат, тригонометрический интеграл и сдвинутое произведение в трёхмерном случае.

Для SEQ прошли все 7 функциональных тестов при запуске:

```bash
PPC_NUM_THREADS=4 OMP_NUM_THREADS=4 ./build/bin/ppc_func_tests \
  --gtest_filter='RectangleIntegrationTests/NazarovaKCalcIntegRectanglesRunFuncTests.*'
```

## 6. Экспериментальная среда

Окружение: Linux 6.6.114.1-microsoft-standard-WSL2, Intel Core i5-1235U, 12 логических CPU, 7.6 GiB RAM, GCC 13.3.0, CMake 3.28.3. Сборка выполнялась в каталоге `build`.

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Releas

cmake --build build -j --parallel

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mpirun --oversubscribe -n 1 ./build/bin/ppc_perf_tests --gtest_filter='*nazarova_k_calc_integ_rectangles_seq*'
```

## 7. Результаты


| Mode     | Workers | Time, s      | Role     |
| -------- | ------- | ------------ | -------- |
| pipeline | 1       | 0.7317971416 | baseline |
| task_run | 1       | 0.6949633072 | baseline |


## 8. Выводы

SEQ-версия задаёт корректный и воспроизводимый baseline. На выбранном performance-тесте она также оказалась самой быстрой, что связано с малой стоимостью вычисления подынтегральной функции по сравнению с накладными расходами параллельных backend-ов.