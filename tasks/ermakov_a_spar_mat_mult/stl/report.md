# Умножение разреженных матриц. Элементы комплексного типа. Формат хранения матрицы – строковый (CRS)

- Студент: Ермаков Алексей Викторович
- Группа: 3823Б1ПР3
- Технология: STL
- Вариант: 6

---

## 1. Назначение реализации

STL-версия реализует многопоточное умножение разреженных матриц
с помощью `std::thread`. В ней разбиение работы между потоками
определяется вручную, без OpenMP и без TBB.

---

## 2. Проверка входных данных

Проверка корректности входных матриц реализована функцией `ValidateMatrix`.
В `ValidationImpl` также проверяется условие совместимости размеров:
`a.cols == b.rows`.

---

## 3. Выбор числа потоков

В этой версии отдельная роль отведена функции `ResolveThreadCount`.
Она принимает:

- число строк результата;
- общую оценку объема работы `total_work`.

Если строк мало или общий объем вычислений меньше `kMinWorkPerThread`,
используется один поток. Иначе число потоков ограничивается сразу тремя факторами:

- значением `ppc::util::GetNumThreads()`;
- числом строк;
- оценкой работы `total_work / kMinWorkPerThread`.

---

## 4. Балансировка нагрузки

Для более равномерного распределения строк между потоками используются
две вспомогательные функции:

- `BuildRowCosts` - оценивает стоимость каждой строки;
- `BuildThreadBoundaries` - по этим стоимостям строит границы диапазонов строк.

Стоимость строки считается как сумма длин строк матрицы `B`,
на которые ссылаются ненулевые элементы строки матрицы `A`.
Таким образом, разбиение учитывает не только число строк,
но и предполагаемую трудоемкость их обработки.

---

## 5. Локальные рабочие структуры

Для потока используется структура `Workspace`, содержащая:

- `accum` - массив накопленных комплексных значений;
- `marks` - массив меток использованных столбцов;
- `touched_cols` - список затронутых столбцов.

Функция `ResetWorkspace` подготавливает рабочее пространство
под текущее число столбцов результата.

Структура `RowData` хранит уже готовые данные строки:

- `cols` - индексы столбцов;
- `vals` - значения элементов строки.

---

## 6. Обработка строки и сборка результата

Функция `MultiplyRow` выполняет полное вычисление одной строки:

- очищает `touched_cols`;
- просматривает ненулевые элементы строки `A`;
- переходит по строкам матрицы `B`;
- накапливает произведения в `accum`;
- сортирует `touched_cols`;
- переносит точные ненулевые значения в `RowData`.

После завершения всех потоков функция `FinalizeResult`:

- строит `row_ptr`;
- вычисляет итоговое число ненулевых элементов;
- копирует `cols` и `vals` всех строк в `c_.col_index` и `c_.values`.

---

## 7. Особенности реализации

- Используется ручной запуск и синхронизация потоков через `std::thread`.
- Разбиение строк учитывает оценку трудоемкости, а не только их количество.
- Для каждой строки применяется та же маркерная схема накопления,
  что и в остальных версиях.
- При `thread_count == 1` версия автоматически переходит к последовательному
  выполнению без создания потоков.

---

## 8. Вывод

STL-реализация отличается от OMP и TBB более явным управлением потоками
и собственной схемой балансировки нагрузки. Именно здесь наиболее явно
проявляется ручное управление распределением строк между потоками.

---

## Приложение (фрагмент кода)

```cpp
void ErmakovASparMatMultSTL::MultiplyRow(int row_index, Workspace &workspace,
                                         RowData &row_data) const {
  workspace.touched_cols.clear();

  for (int ak = a_.row_ptr[row_index]; ak < a_.row_ptr[row_index + 1]; ++ak) {
    const int b_row = a_.col_index[ak];
    const auto a_value = a_.values[ak];

    for (int bk = b_.row_ptr[b_row]; bk < b_.row_ptr[b_row + 1]; ++bk) {
      const int col = b_.col_index[bk];
      const auto product = a_value * b_.values[bk];

      if (workspace.marks[static_cast<std::size_t>(col)] != row_index) {
        workspace.marks[static_cast<std::size_t>(col)] = row_index;
        workspace.accum[static_cast<std::size_t>(col)] = product;
        workspace.touched_cols.push_back(col);
      } else {
        workspace.accum[static_cast<std::size_t>(col)] += product;
      }
    }
  }

  std::ranges::sort(workspace.touched_cols);

  row_data.cols.clear();
  row_data.vals.clear();
  row_data.cols.reserve(workspace.touched_cols.size());
  row_data.vals.reserve(workspace.touched_cols.size());

  for (int col : workspace.touched_cols) {
    const auto &value = workspace.accum[static_cast<std::size_t>(col)];
    if (value != kZero) {
      row_data.cols.push_back(col);
      row_data.vals.push_back(value);
    }
  }
}
```
