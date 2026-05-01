# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – строковый (CRS)

- Студент: Ермаков Алексей Викторович
- Группа: 3823Б1ПР3
- Технология: OMP
- Вариант: 6

---

## 1. Назначение реализации

OpenMP-версия распараллеливает вычисление строк результирующей матрицы между
потоками, сохраняя ту же математическую схему, что и последовательная версия.

---

## 2. Проверка входных данных

Проверка входных матриц полностью повторяет последовательную версию.
Используется функция `ValidateMatrix`, а в `ValidationImpl` дополнительно
контролируется условие `a.cols == b.rows`.

---

## 3. Схема распараллеливания

Распараллеливание организовано по строкам матрицы `A`.

В `RunImpl` создаются два массива построчных результатов:

- `row_values` - ненулевые значения для каждой строки;
- `row_cols` - соответствующие индексы столбцов.

Далее выполняется `#pragma omp parallel`, внутри которого каждый поток создает
свои локальные рабочие структуры:

- `row_vals` - массив частичных сумм по столбцам;
- `row_mark` - массив меток использованных столбцов;
- `used_cols` - список затронутых столбцов.

Цикл по строкам распараллеливается директивой `#pragma omp for`.
Каждая строка обрабатывается независимо.

---

## 4. Локальная обработка строки

Для вычисления строки используются вспомогательные функции:

- `AccumulateRowProducts` - накапливает произведения `A[i,j] * B[j,k]`;
- `SortUsedCols` - сортирует список использованных столбцов;
- `CollectRowValues` - переносит точные ненулевые значения строки
  в контейнеры `cols` и `vals`.

В отличие от последовательной версии, здесь результат строки не записывается
сразу в `c_`, а сначала сохраняется в `row_cols[i]` и `row_values[i]`.

---

## 5. Сборка результата

После завершения параллельной секции выполняется последовательная сборка матрицы `c_`:

1. Проходом по `row_values` строится массив `row_ptr`.
2. Подсчитывается общее число ненулевых элементов `nnz`.
3. `values` и `col_index` резервируют нужную память.
4. Данные строк по порядку вставляются в `c_.values` и `c_.col_index`.

---

## 6. Особенности реализации

- Каждый поток работает только со своими локальными векторами.
- Фильтрация нулей выполняется в `CollectRowValues` по точному сравнению
  с комплексным нулем.
- Сборка итогового CRS вынесена из параллельной части, что упрощает
  формирование корректного массива `row_ptr`.

---

## 7. Вывод

OMP-реализация представляет собой прямое распараллеливание базового алгоритма
по строкам. Она сохраняет ту же логику вычислений, но переносит обработку строк
в независимые потоки OpenMP и затем последовательно собирает итоговый результат.

---

## Приложение (фрагмент кода)

```cpp
#pragma omp parallel default(none) shared(m, p, row_values, row_cols)
{
  std::vector<std::complex<double>> row_vals(
      static_cast<std::size_t>(p), std::complex<double>(0.0, 0.0));
  std::vector<int> row_mark(static_cast<std::size_t>(p), -1);
  std::vector<int> used_cols;
  used_cols.reserve(256);

#pragma omp for
  for (int i = 0; i < m; ++i) {
    AccumulateRowProducts(i, row_vals, row_mark, used_cols);
    SortUsedCols(used_cols);

    const auto row_i = static_cast<std::size_t>(i);

    CollectRowValues(row_vals, used_cols, row_cols[row_i], row_values[row_i]);
  }
}

int nnz = 0;

for (int i = 0; i < m; ++i) {
  const auto row_i = static_cast<std::size_t>(i);
  c_.row_ptr[row_i] = nnz;
  nnz += static_cast<int>(row_values[row_i].size());
}

c_.row_ptr[static_cast<std::size_t>(m)] = nnz;
c_.values.reserve(static_cast<std::size_t>(nnz));
c_.col_index.reserve(static_cast<std::size_t>(nnz));

for (int i = 0; i < m; ++i) {
  const auto row_i = static_cast<std::size_t>(i);
  CollectRowValues(row_vals, used_cols, row_cols[static_cast<std::size_t>(i)],
                   row_values[static_cast<std::size_t>(i)]);
}
```
