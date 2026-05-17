# Умножение матриц: алгоритм Кэннона — OMP

- Student: Зязева Светлана Александровна
- Technology: OMP
- Variant: zyazeva_s_matrix_mult_cannon_alg

## 1. Контекст

OpenMP-версия реализует алгоритм Кэннона для блочного умножения матриц,
распараллеливая независимые блочные операции между потоками. В отличие от
SEQ, где используется простой тройной цикл, здесь матрицы разбиваются на
блоки по сетке `grid_size × grid_size`, и шаги алгоритма Кэннона
выполняются параллельно. При невозможности применить Кэннона (число
потоков не является полным квадратом, или размер матрицы не кратен
`sqrt(num_threads)`) реализация откатывается к параллельному тройному
циклу.

## 2. Постановка задачи

**Входные данные**
(`InType = std::tuple<size_t, std::vector<double>, std::vector<double>>`):

- `sz` — сторона квадратной матрицы (тип `size_t`);
- `m1` — матрица A, построчно, размер `sz * sz`;
- `m2` — матрица B, аналогично.

**Выходные данные** (`OutType = std::vector<double>`):
результирующая матрица C = A·B, построчно, размер `sz * sz`.

**Ограничения:**

- `sz > 0` — нулевой размер недопустим;
- `m1.size() == sz * sz` и `m2.size() == sz * sz` — проверяется в
  `ValidationImpl`;
- элементы — `double`, переполнение не рассматривается.

**Крайние случаи:**

- матрица 1×1 — корректно обрабатывается основным циклом;
- пустые векторы при `sz > 0` — отклоняются валидацией.

## 3. Базовый алгоритм

Алгоритм Кэннона разбивает матрицы A и B на `g × g` блоков
(где `g = sqrt(num_threads)`). Начальное смещение расставляет блоки A
по строкам, а B по столбцам. Затем выполняется `g` шагов: на каждом шаге
блоки перемножаются и сдвигаются (A — влево, B — вверх). Результирующие
блоки собираются в итоговую матрицу. Метод `RegularMultiplication`
применяется как fallback.

## 4. Схема распараллеливания

**Параллелизуемые области:**

1. `AlignBlocks` — начальное выравнивание блоков;
2. `CannonStep` — умножение и сдвиг блоков на каждом шаге;
3. `AssembleResult` — сборка итоговой матрицы из блоков;
4. `RegularMultiplication` (fallback) — цикл по строкам.

```cpp
// File: omp/src/ops_omp.cpp — AlignBlocks
#pragma omp parallel for default(none) \
    shared(blocks_a, blocks_b, aligned_a, aligned_b, \
           grid_size, grid_size_t) \
    collapse(2)
for (int i = 0; i < grid_size; ++i) {
  for (int j = 0; j < grid_size; ++j) {
    // aligned_a[block_idx] = blocks_a[...];
    // aligned_b[block_idx] = blocks_b[...];
  }
}
```

**Атрибуты переменных:**

- `default(none)` — все атрибуты перечислены явно;
- `shared`: блочные векторы (`blocks_a`, `blocks_b`, `aligned_a`,
  `aligned_b`, `blocks_c`, `res_m`), размерные константы (`grid_size`,
  `block_size` и их `size_t`-аналоги);
- `private`: итераторы `i` и `j` неявно private при `parallel for`;

**`collapse(2)`** применён к двойному циклу по `grid_size × grid_size`.
Это позволяет OpenMP распределять все `grid_size²` итераций, а не только
строки внешнего цикла — важно при малом `grid_size`.

**Барьеры:** в конце каждой `parallel for` стоит неявный барьер. Между
шагами Кэннона барьер выполняется при завершении `CannonStep`: следующий
шаг не начнётся раньше, чем все блоки перемножены и сдвинуты.

**Число потоков:** задаётся переменной `OMP_NUM_THREADS`
(= `PPC_NUM_THREADS`).

## 5. Детали реализации

**Файлы:** `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`

`RunImpl` определяет применимость алгоритма Кэннона проверкой
`IsPerfectSquare(num_threads) && sz >= num_threads`
`&& (sz % (int)sqrt(num_threads) == 0)`.

```cpp
// File: omp/src/ops_omp.cpp — CannonStep
#pragma omp parallel for default(none) \
    shared(aligned_a, aligned_b, blocks_c, \
           grid_size, block_size, grid_size_t) \
    collapse(2)
for (int i = 0; i < grid_size; ++i) {
  for (int j = 0; j < grid_size; ++j) {
    const size_t block_idx =
        (static_cast<size_t>(i) * grid_size_t)
        + static_cast<size_t>(j);
    MultiplyBlocks(
        aligned_a[block_idx], aligned_b[block_idx],
        blocks_c[block_idx], block_size);
  }
}
```

Каждый поток работает со своим `block_idx`, записывая в
`blocks_c[block_idx]`. Блоки не пересекаются — `atomic` и `critical`
не нужны. Сдвиг выполняется параллельно: `new_aligned_*` создаются
локально и заменяют старые после неявного барьера.

Риски гонок устранены: `MultiplyBlocks` обновляет только
`blocks_c[block_idx]`.

## 6. Проверка корректности

OMP-версия сравнивалась с SEQ для всех функциональных тестов, включая
матрицы 1×1, 2×2 и средних размеров. Ветка `RegularMultiplication`
тестировалась при `PPC_NUM_THREADS = 1` и нечётных значениях.
Расхождений не наблюдалось. Дополнительно проверялась устойчивость при
`PPC_NUM_THREADS = 1, 2, 4`.

## 7. Экспериментальная среда

| Параметр | Значение |
| --- | --- |
| CPU | AMD Ryzen 7 8845HS w/ Radeon 780M Graphics |
| RAM | 32 ГБ |
| OS | Windows 11 Pro |
| Compiler | g++ / clang++ (Release) |
| CMake build type | Release |
| OMP_NUM_THREADS | задаётся runner-ом курса |

**Команды запуска:**

```bash
    export OMP_NUM_THREADS=4

    ./build/bin/ppc_func_tests --gtest_filter=*"zyazeva_s"*

    ./build/bin/ppc_perf_tests --gtest_filter=*"zyazeva_s"*
```

## 8. Результаты

| Режим | Время OMP (с) | Время SEQ (с) | Ускорение |
| --- | --- | --- | --- |
| pipeline | 0.1743746712 | 0.7128969324 | **4.09x** |
| task_run | 0.2143771898 | 0.5624887756 | **2.62x** |

Ускорение в режиме `pipeline` выше, чем в `task_run`. В обоих случаях
OMP даёт значительный прирост. При увеличении числа потоков сверх
`grid_size^2` ускорение ограничено алгоритмически.

## 9. Выводы

OpenMP дал ускорение 2.6–4.1x. Ключевой фактор — параллельная обработка
независимых блоков алгоритма Кэннона. Ограничения: число потоков должно
быть полным квадратом и делителем `sz`; иначе код откатывается к
`RegularMultiplication`. Overhead на аллокацию `new_aligned_*` на каждом
шаге ограничивает дальнейшее масштабирование.
