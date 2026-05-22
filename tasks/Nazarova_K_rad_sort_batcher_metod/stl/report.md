# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) — STL

- Student: Назарова Ксения Олеговна
- Technology: STL
- Variant: 9

## 1. Контекст

STL-версия показывает ручное разбиение линейного диапазона ячеек между асинхронными задачами стандартной библиотеки C++. В отличие от OpenMP и TBB, здесь границы поддиапазонов вычисляются явно в коде.

## 2. Постановка задачи

Нужно вычислить тот же многомерный интеграл методом средних прямоугольников, что и в SEQ-версии. Последовательная реализация используется как эталон результата и baseline для сравнения времени.

## 3. Базовый алгоритм

Сначала вычисляется общее число ячеек `total_cells_`. Далее диапазон `[0, total_cells_)` делится на почти равные чанки. Каждая асинхронная задача обходит свой чанк, восстанавливает координаты центра ячейки и возвращает локальную сумму.

## 4. Схема распараллеливания

Число задач берётся из `std::thread::hardware_concurrency()`, затем ограничивается числом ячеек. Размер базового чанка равен `total_cells_ / thread_count`, остаток распределяется по первым задачам.

Ключевой фрагмент из `stl/src/ops_stl.cpp`:

```cpp
std::vector<std::future<double>> futures;
futures.reserve(thread_count);

std::size_t begin = 0U;
for (std::size_t thread_id = 0U; thread_id < thread_count; ++thread_id) {
  const std::size_t extra = thread_id < remainder ? 1U : 0U;
  const std::size_t end = begin + base_chunk + extra;
  futures.emplace_back(std::async(std::launch::async, [&, begin, end] {
    std::vector<double> point(dimension_, 0.0);
    double local_sum = 0.0;
    for (std::size_t linear_index = begin; linear_index < end; ++linear_index) {
      FillPointFromLinearIndex(input, step_sizes_, linear_index, point);
      local_sum += input.function(point);
    }
    return local_sum;
  }));
  begin = end;
}
```

Синхронизация общей суммы через mutex не нужна: каждая задача возвращает своё значение, а итоговая сумма собирается в основном потоке после `future.get()`. Локальный вектор `point` также принадлежит только одной задаче.

## 5. Детали реализации

Файлы реализации: `stl/include/ops_stl.hpp`, `stl/src/ops_stl.cpp`.

`ValidationImpl` проверяет корректность входа и переполнение `total_cells`. `PreProcessingImpl` вычисляет размеры шагов, объём ячейки и общее число линейных итераций. В `RunImpl` создаются асинхронные задачи, после чего основной поток последовательно получает результаты через `future.get()` и суммирует их.

## 6. Проверка корректности

STL-версия прошла все 7 функциональных тестов общего набора. Поскольку общих изменяемых буферов между задачами нет, гонка возможна только при небезопасной пользовательской функции; тестовые функции чистые и не используют разделяемое состояние.

## 7. Экспериментальная среда

Окружение: Linux 6.6.114.1-microsoft-standard-WSL2, Intel Core i5-1235U, 12 логических CPU, GCC 13.3.0. Основной запуск выполнялся с `PPC_NUM_THREADS=4`, но сама STL-реализация ориентируется на `std::thread::hardware_concurrency()`.

Команды запуска:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Releas

cmake --build build -j --parallel

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

mpirun --oversubscribe -n 4 ./build/bin/ppc_perf_tests --gtest_filter='*nazarova_k_calc_integ_rectangles_stl*'
```

## 8. Результаты

| Mode | Workers in report | Time, s | Speedup | Efficiency |
|---|---:|---:|---:|---:|
| pipeline | 4 | 2.0161764572 | 0.36 | 0.09 |
| task_run | 4 | 1.9265005476 | 0.36 | 0.09 |

STL-версия оказалась быстрее OMP и TBB в этом запуске, но всё равно медленнее SEQ. Главный источник потерь — запуск асинхронных задач и сбор `future`, которые дороги относительно простого тела цикла.

## 9. Выводы

Ручное разбиение диапазона реализовано корректно: поддиапазоны не пересекаются, локальные суммы независимы, объединение выполняется после завершения задач. Для более тяжёлого интегранта или при явном управлении числом задач эта схема может быть эффективнее, но текущий тест не компенсирует накладные расходы `std::async`.
