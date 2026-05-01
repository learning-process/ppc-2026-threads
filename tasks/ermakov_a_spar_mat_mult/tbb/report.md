# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – строковый (CRS)

- Студент: Ермаков Алексей Викторович
- Группа: 3823Б1ПР3
- Технология: TBB
- Вариант: 6

---

## 1. Назначение реализации

TBB-версия использует `oneapi::tbb::parallel_for` для параллельной обработки
диапазонов строк матрицы `A`. По математической логике она совпадает с SEQ и OMP,
но отличается способом организации локальной памяти и задач.

---

## 2. Проверка входных данных

Проверка входных матриц реализована так же, как и в SEQ/OMP:

- используется функция `ValidateMatrix`;
- в `ValidationImpl` проверяется условие `a.cols == b.rows`.

---

## 3. Организация локальной памяти

В этой версии введена локальная структура `RowWorkspace`, которая содержит:

- `row_vals` - массив частичных сумм по столбцам;
- `row_mark` - массив меток использованных столбцов;
- `used_cols` - список затронутых столбцов.

Экземпляры `RowWorkspace` создаются через
`tbb::enumerable_thread_specific<RowWorkspace>`, поэтому каждый рабочий поток TBB
получает собственный набор буферов.

---

## 4. Схема распараллеливания

В `RunImpl` строки матрицы `A` делятся на блоки с помощью
`tbb::blocked_range<int>(0, m, grain_size)`.

Размер блока выбирается функцией `ResolveGrainSize`, которая вычисляет
приблизительно `rows / 16`, но не меньше единицы. Идея в том, чтобы получить
достаточное число задач без излишнего дробления.

Внутри `parallel_for` для каждой строки вызываются:

- `AccumulateRowProducts` - накопление произведений;
- `SortUsedCols` - сортировка использованных столбцов;
- `CollectRowValues` - перенос точных ненулевых элементов в построчные контейнеры.

---

## 5. Сборка результата

После завершения `parallel_for` выполняется последовательная сборка матрицы `c_`:

- проходом по `row_values` вычисляется `row_ptr`;
- резервируется память под итоговые `values` и `col_index`;
- строки по очереди копируются в итоговый CRS.

Таким образом, TBB отвечает только за вычисление построчных данных,
а финальная компоновка структуры CRS остается последовательной.

---

## 6. Особенности реализации

- Локальные рабочие буферы создаются через `enumerable_thread_specific`,
  а не вручную.
- Используется явный `grain_size`, зависящий от числа строк.
- Логика фильтрации нулей и порядок обработки строки совпадают
  с остальными shared-memory версиями.

---

## 7. Вывод

TBB-реализация использует задачно-ориентированную схему параллелизма:
строки матрицы делятся на диапазоны, а локальная память хранится
в поток-специфичных рабочих структурах. Это делает реализацию компактной
и хорошо согласованной с моделью TBB.

---

## Приложение (фрагмент кода)

```cpp
struct RowWorkspace {
  std::vector<std::complex<double>> row_vals;
  std::vector<int> row_mark;
  std::vector<int> used_cols;

  explicit RowWorkspace(int cols)
      : row_vals(static_cast<std::size_t>(cols), std::complex<double>(0.0, 0.0)),
        row_mark(static_cast<std::size_t>(cols), -1) {
    used_cols.reserve(256);
  }
};

int ResolveGrainSize(int rows) {
  if (rows <= 0) {
    return 1;
  }

  constexpr int kTargetChunks = 16;
  return std::max(1, rows / kTargetChunks);
}

tbb::enumerable_thread_specific<RowWorkspace> workspace([&] { return RowWorkspace(p); });
const int grain_size = ResolveGrainSize(m);

tbb::parallel_for(tbb::blocked_range<int>(0, m, grain_size), [&](const tbb::blocked_range<int> &range) {
  auto &local = workspace.local();

  for (int i = range.begin(); i != range.end(); ++i) {
    AccumulateRowProducts(i, local.row_vals, local.row_mark, local.used_cols);
    SortUsedCols(local.used_cols);

    const auto row_i = static_cast<std::size_t>(i);
    CollectRowValues(local.row_vals, local.used_cols, row_cols[row_i], row_values[row_i]);
  }
});
```
