# Вычисление многомерных интегралов с использованием многошаговой схемы (метод прямоугольников) — ALL

- Student: Назарова Ксения Олеговна
- Technology: ALL, MPI + TBB
- Variant: 9

## 1. Контекст

Гибридная версия сочетает межпроцессное разбиение диапазона через MPI и внутрипроцессную параллельную редукцию через
oneTBB. Такая схема проверяет два уровня параллелизма: распределение ячеек между rank-ами и обработку локального
диапазона внутри каждого процесса.

## 2. Постановка задачи

Постановка совпадает с общей задачей: требуется вычислить многомерный интеграл методом средних прямоугольников.
SEQ-версия является эталоном результата, а корневой отчёт содержит сводное сравнение всех backend-ов.

## 3. Базовый алгоритм

Все ячейки сетки нумеруются линейно от `0` до `total_cells_ - 1`. MPI-процессы получают непересекающиеся отрезки
этого диапазона. Каждый процесс считает локальную сумму на своём отрезке, после чего глобальная сумма собирается
коллективной операцией.

## 4. Межпроцессная схема

В начале `RunImpl` процесс узнаёт свой ранг и размер коммуникатора:

```cpp
int rank = 0;
int size = 1;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
```

Далее диапазон делится почти поровну. Остаток от деления распределяется между первыми rank-ами, поэтому размеры
локальных чанков отличаются не более чем на одну ячейку. После локального вычисления используется `MPI_Allreduce`,
чтобы каждый процесс получил один и тот же итоговый результат.

```cpp
double global_sum = 0.0;
MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

result_ = global_sum * cell_volume_;
GetOutput() = result_;
```

`MPI_Allreduce` одновременно выполняет суммирование и распространение результата всем участникам коммуникатора. Это
важно для тестовой инфраструктуры, где результат может проверяться на каждом rank-е.

## 5. Внутрипроцессная схема

Внутри локального диапазона используется `tbb::parallel_reduce`:

```cpp
const double local_sum = tbb::parallel_reduce(tbb::blocked_range<std::size_t>(local_begin, local_end), 0.0,
                                              [&](const tbb::blocked_range<std::size_t> &range, double partial_sum) {
  std::vector<double> point(dimension_, 0.0);
  for (std::size_t linear_index = range.begin(); linear_index < range.end(); ++linear_index) {
    FillPointFromLinearIndex(input, step_sizes_, linear_index, point);
    partial_sum += input.function(point);
  }
  return partial_sum;
}, std::plus<>());
```

Каждый TBB-поддиапазон имеет собственный вектор координат и собственную частичную сумму. Совместная запись отсутствует
до MPI-редукции.

## 6. Детали реализации

Файлы реализации: `all/include/ops_all.hpp`, `all/src/ops_all.cpp`.

`ValidationImpl` проверяет вход и переполнение общего числа ячеек. `PreProcessingImpl` вычисляет шаги, объём ячейки и
`total_cells_`. В `RunImpl` происходит двухуровневое разбиение: сначала MPI делит линейный диапазон между процессами,
затем TBB делит локальный диапазон внутри процесса.

Потенциальные узкие места: коммуникация `MPI_Allreduce`, oversubscription при слишком большом произведении
`ranks × threads`, а также overhead TBB для лёгкой функции.

## 7. Проверка корректности

Обычный запуск без `mpirun` пропускает `all`-тесты, потому что MPI-задачи должны исполняться под MPI launcher.
Отдельный запуск:

```bash
mpirun --allow-run-as-root -np 2 ./build/bin/ppc_func_tests \
  --gtest_filter='RectangleIntegrationTests/NazarovaKCalcIntegRectanglesRunFuncTests.*all_enabled*'
```

Результат: все 7 функциональных тестов `all` прошли. Это подтверждает, что распределение диапазона между двумя
rank-ами и последующая `MPI_Allreduce` дают тот же результат, что и эталонные значения.

## 8. Экспериментальная среда

Окружение: Linux 6.6.114.1-microsoft-standard-WSL2, Intel Core i5-1235U, 12 логических CPU, GCC 13.3.0. Для гибридного
запуска использовалась конфигурация `2 × 4`: два MPI-процесса и `PPC_NUM_THREADS=4`.

Команда запуска для гибридной конфигурации:

```bash
cmake -S . -B build -D USE_COVERAGE=ON -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D CMAKE_BUILD_TYPE=Release

cmake --build build -j --parallel

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

PPC_NUM_PROC=2 PPC_NUM_THREADS=4 \
mpirun --oversubscribe -x PPC_NUM_PROC -x PPC_NUM_THREADS -n 2 \
  ./build/bin/ppc_perf_tests --gtest_filter='*nazarova_k_calc_integ_rectangles_all*'
```

## 9. Результаты

| Mode     | Ranks | Threads per rank | Total workers | Time, s      | Speedup | Efficiency |
| -------- | ----- | ---------------- | ------------- | ------------ | ------- | ---------- |
| pipeline | 2     | 4                | 8             | 1.0314585568 | 0.71    | 0.09       |
| task_run | 2     | 4                | 8             | 1.0240943998 | 0.68    | 0.08       |

Для сравнения, однопроцессный `all` в общем performance-запуске показал 1.9185 s для `pipeline` и 1.9017 s для
`task_run`. Переход к двум MPI-процессам уменьшил время, но не дал ускорения относительно SEQ.

## 10. Выводы

ALL-версия корректно реализует иерархическое распараллеливание: MPI делит диапазон между процессами, TBB считает
локальные суммы, `MPI_Allreduce` объединяет результат. На выбранном тесте гибридная схема снижает время относительно
однопроцессного `all`, но эффективность остаётся низкой из-за стоимости TBB runtime, MPI-редукции и лёгкого тела
интегрируемой функции.
