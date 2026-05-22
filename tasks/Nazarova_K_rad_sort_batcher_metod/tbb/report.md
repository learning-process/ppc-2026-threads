# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) — TBB

- Student: Назарова Ксения Олеговна
- Technology: TBB
- Variant: 9

## 1. Контекст

TBB-версия использует задачно-ориентированную модель параллелизма. Поскольку итог интегрирования является суммой независимых вкладов ячеек, естественным примитивом является `tbb::parallel_reduce`.

## 2. Постановка задачи

Задача совпадает с последовательной постановкой: вычислить приближённый многомерный интеграл методом средних прямоугольников. SEQ-версия служит эталоном результата и baseline для производительности.

## 3. Базовый алгоритм

Как и в OMP-версии, многомерная сетка линеаризуется в диапазон `[0, total_cells_)`. Для каждого линейного индекса восстанавливаются координатные индексы, затем строится центр ячейки и вычисляется значение функции.

## 4. Схема распараллеливания

Основной диапазон задаётся как `tbb::blocked_range<std::size_t>(0U, total_cells_)`. Runtime oneTBB разбивает этот диапазон на поддиапазоны, каждый worker накапливает локальную сумму, а затем частичные результаты объединяются через `std::plus<>`.

Ключевой фрагмент из `tbb/src/ops_tbb.cpp`:

```cpp
const double sum = tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0U, total_cells_), 0.0,
                                        [&](const tbb::blocked_range<std::size_t> &range, double local_sum) {
  std::vector<double> point(dimension_, 0.0);
  for (std::size_t linear_index = range.begin(); linear_index < range.end(); ++linear_index) {
    FillPointFromLinearIndex(input, step_sizes_, linear_index, point);
    local_sum += input.function(point);
  }
  return local_sum;
}, std::plus<>());
```

Явный `grainsize` и partitioner не задаются, поэтому используется стандартная стратегия разбиения oneTBB. Конкуренция в локальном запуске ограничивалась окружением тестового runner-а через `PPC_NUM_THREADS=4`.

## 5. Детали реализации

Файлы реализации: `tbb/include/ops_tbb.hpp`, `tbb/src/ops_tbb.cpp`.

Функция `FillPointFromLinearIndex` вынесена отдельно, чтобы одинаково преобразовывать линейный индекс в координаты центра ячейки. Внутри тела `parallel_reduce` вектор `point` создаётся локально для каждого поддиапазона, поэтому совместной записи в него нет. Общая сумма не изменяется напрямую: каждый task возвращает локальный результат, а TBB объединяет его редукцией.

## 6. Проверка корректности

TBB-версия прошла все 7 функциональных тестов из общего набора. Проверялись одномерные, двумерные и трёхмерные случаи, включая тригонометрический интеграл с ненулевой численной погрешностью.

## 7. Экспериментальная среда

Окружение: Linux 6.6.114.1-microsoft-standard-WSL2, Intel Core i5-1235U, 12 логических CPU, GCC 13.3.0. Основная конфигурация: `PPC_NUM_THREADS=4`.

Команды запуска:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Releas

cmake --build build -j --parallel

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mpirun --oversubscribe -n 4 ./build/bin/ppc_perf_tests --gtest_filter='*nazarova_k_calc_integ_rectangles_tbb*'
```

## 8. Результаты

| Mode | Workers | Time, s | Speedup | Efficiency |
|---|---:|---:|---:|---:|
| pipeline | 4 | 1.1344363388 | 0.65 | 0.16 |
| task_run | 4 | 1.0996614686 | 0.63 | 0.16 |

На выбранном тесте TBB не даёт ускорения относительно SEQ. Основной фактор — overhead `parallel_reduce` и разбиения диапазонов при дешёвом вычислении функции.

## 9. Выводы

TBB-реализация корректно выражает задачу как параллельную редукцию по независимому диапазону ячеек. Подход удобен и безопасен по данным, но для лёгкого интегранта и текущего размера сетки runtime-издержки оказываются больше полезного выигрыша.
