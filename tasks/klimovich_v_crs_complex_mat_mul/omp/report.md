# Умножение разреженных матриц комплексного типа в формате CRS — OMP

- Student: Климович Виктор Олегович, group 3823Б1ПР4
- Technology: OpenMP
- Variant: 6

## 1. Контекст

OMP-версия переносит row-wise Густавсона из SEQ в параллельное окружение
OpenMP. Поскольку обработка разных строк `i` матрицы `A` независима (каждая
формирует свою строку `C[i, :]`), это естественный data-parallel случай.
Цель отчёта — показать, что независимость данных не нарушена и при этом
каждый поток имеет приватный SPA достаточного размера, чтобы исключить гонки
и false sharing на ключевом аккумуляторе.

## 2. Постановка задачи

См. [seq/report.md](../seq/report.md). Вход, выход, ограничения и численный
критерий — те же.

## 3. Базовый алгоритм

Тот же Густавсон row-wise со SPA. Изменяется только распределение работы:
итерации по `i` (внешнему циклу строк) распределяются между потоками
OpenMP, каждый поток держит **свой собственный** `spa`, `touched_by_row`,
`touched_cols`.

## 4. Схема распараллеливания

- Параллелится **внешний цикл по строкам**
  `for (int i = 0; i < lhs.n_rows; ++i)`.
- На уровне `#pragma omp parallel` объявлены приватные буферы внутри блока —
  они автоматически private для каждого потока, потому что объявлены в
  области видимости потока.
- `lhs`, `rhs`, `per_row` — `shared`. `lhs` и `rhs` только читаются;
  `per_row[i]` пишется каждым потоком, но **каждый поток пишет в свою
  собственную ячейку** (свою строку), поэтому записи независимы.
- `schedule(dynamic, 16)` — динамическое распределение чанками по 16 строк.
  Выбран `dynamic`, а не `static`, потому что число ненулей в строке
  (а значит и стоимость её обработки) сильно колеблется между строками;
  чанк 16 — компромисс между накладными расходами на выдачу заданий и
  балансом нагрузки.
- Аккумуляция результата выполнена через **двухфазную сборку**: сначала
  каждый поток складывает результат в `per_row[i]` (плотные временные буферы
  для каждой строки), затем вне параллельной области функция `Assemble`
  собирает итоговый CRS. Это позволяет полностью избежать критических
  секций.
- В конце `#pragma omp parallel` есть неявный барьер — он необходим, чтобы
  `Assemble` стартовал только после того, как все строки готовы.

## 5. Детали реализации

Файлы: [omp/include/ops_omp.hpp](include/ops_omp.hpp),
[omp/src/ops_omp.cpp](src/ops_omp.cpp).

Параллельная область:

```cpp
// File: omp/src/ops_omp.cpp
#pragma omp parallel default(none) shared(lhs, rhs, per_row)
{
  std::vector<Cplx> spa(static_cast<std::size_t>(rhs.n_cols));
  std::vector<int> touched_by_row(static_cast<std::size_t>(rhs.n_cols), -1);
  std::vector<int> touched_cols;
  touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));

#pragma omp for schedule(dynamic, 16)
  for (int i = 0; i < lhs.n_rows; ++i) {
    GustavsonRow(lhs, rhs, i, spa, touched_by_row, touched_cols, per_row[i]);
  }
}
```

`default(none)` обязателен по правилам clang-tidy
(`openmp-use-default-none`) — заставляет явно перечислить, какие переменные
shared. Все буферы, объявленные **внутри** `parallel`-блока, автоматически
приватны: это и есть наш SPA, sentinel и список затронутых столбцов.
Поэтому никакой `private(...)` оговорки не требуется.

Важное замечание про MSVC C3052: значение `Cplx(0.0, 0.0)` как явный вызов
конструктора внутри `default(none)`-области ломает MSVC OpenMP (тип `Cplx`
рассматривается как имя, требующее data-sharing). Поэтому здесь
используется `std::vector<Cplx>` без второго аргумента —
value-initialization элементов комплексного типа даёт ровно `(0, 0)`, без
явной ссылки на `Cplx` как на выражение.

Финальная сборка после барьера:

```cpp
// File: omp/src/ops_omp.cpp
CrsMatrix Assemble(int rows, int cols, const std::vector<RowStage> &per_row) {
  CrsMatrix out(rows, cols);
  for (int i = 0; i < rows; ++i) {
    out.row_offsets[i + 1] = out.row_offsets[i]
        + static_cast<int>(per_row[i].cols.size());
  }
  out.col_indices.reserve(static_cast<std::size_t>(out.row_offsets[rows]));
  out.data.reserve(static_cast<std::size_t>(out.row_offsets[rows]));
  for (int i = 0; i < rows; ++i) {
    out.col_indices.insert(out.col_indices.end(),
                           per_row[i].cols.begin(), per_row[i].cols.end());
    out.data.insert(out.data.end(),
                    per_row[i].vals.begin(), per_row[i].vals.end());
  }
  return out;
}
```

Эта стадия последовательна, но дёшева: O(nnz(C)) одной линейной прохода без
случайного доступа. Альтернатива с записью напрямую в общий CRS требовала
бы либо префиксной суммы и двух проходов, либо критических секций — мы
намеренно избегаем и того и другого.

## 6. Проверка корректности

Functional-тест [tests/functional/main.cpp](../tests/functional/main.cpp)
запускает OMP-версию через `AddFuncTask<KlimovichVCrsComplexMatMulOmp, InType>`
на тех же 10 наборах размеров, что и SEQ, со сравнением с плотным эталоном
`DenseMultiply` и допуском `1e-9`. Дополнительно проверено, что результат
численно совпадает с SEQ для всех `PPC_NUM_THREADS ∈ {1, 2, 4, 8}`.

Гонок нет по построению: каждый поток работает с приватным `spa`, и каждая
запись в `per_row[i]` происходит ровно один раз и только в «своей» ячейке
итерации `i`.

## 7. Экспериментальная среда

- CPU: Intel Core i7-11800H @ 2.30 GHz (8C/16T)
- RAM: 16 GB DDR4-3200
- OS: Windows 11 Pro 22H2 (build 22631)
- Compiler: MSVC 19.44.35211 (Visual Studio 2022)
- Build type: Release
- `PPC_NUM_THREADS` ∈ {1, 2, 4, 8}; `OMP_NUM_THREADS` экспортируется
  runner-ом курса автоматически.

Команды:

```bash
set PPC_NUM_THREADS=4
scripts/run_tests.py --running-type=threads --counts 1 2 4 8
scripts/run_tests.py --running-type=performance
```

## 8. Результаты

Размер задачи `1500×1500`, `T_seq(task) = 0.152 s`,
`T_seq(pipeline) = 1.498 s`. Медианы по 10 повторам.

| threads | mode     | median time, s | speedup vs seq | efficiency, % |
| ------: | -------- | -------------: | -------------: | ------------: |
|       1 | task     |          0.155 |           0.98 |            98 |
|       2 | task     |          0.082 |           1.85 |            93 |
|       4 | task     |          0.044 |           3.45 |            86 |
|       8 | task     |          0.029 |           5.24 |            66 |
|       1 | pipeline |          1.522 |           0.98 |            98 |
|       2 | pipeline |          0.803 |           1.87 |            93 |
|       4 | pipeline |          0.418 |           3.58 |            90 |
|       8 | pipeline |          0.272 |           5.51 |            69 |

Ускорение `S(p) = T_seq / T(p)`, эффективность `E(p) = S(p) / p`.

Видимые закономерности:

- До 4 потоков масштабирование почти линейное (efficiency 86–93%): задача
  хорошо разделяется по строкам, накладные расходы OMP-runtime малы.
- При переходе с 4 на 8 потоков эффективность падает с ≈86–90% до 66–69%.
  Причины: (а) на ноутбуке 8 физических ядер делят общую LLC и шину памяти,
  и алгоритм memory-bound на доступах в `spa`; (б) суммарный footprint
  приватных SPA становится `8 * 1500 * 16 B ≈ 188 KB`, что уже сопоставимо
  с типовым размером L2 на ядро.
- Pipeline-режим стабильнее task-режима (выше эффективность при тех же
  потоках) — это ожидаемо: повторные прогоны прогревают кэш и амортизируют
  стартовые издержки.

## 9. Выводы

OMP-версия даёт стабильное ускорение на матрице `1500×1500` за счёт полной
независимости работ по строкам. На 4 потоках достигается ≈3.5× ускорение
при 86–90% эффективности; на 8 потоках суммарное ускорение растёт до ≈5.5×,
но эффективность падает из-за memory-bandwidth bottleneck. Главные
методические решения: двухфазная сборка через `per_row[i]` (исключает
критические секции и синхронизацию на общем CRS) и приватный SPA на поток
(исключает false sharing).
