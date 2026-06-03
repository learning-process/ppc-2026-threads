# Умножение разреженных матриц с комплексными элементами - OMP

- **Student:** Люлин Ярослав Сергеевич
- **Technology:** OpenMP (OMP)
- **Variant:** 7

## 1. Контекст

OpenMP-версия является первым этапом в распараллеливании последовательного алгоритма
умножения матриц в формате CCS. Основная цель — оценить эффективность директивного подхода
для задачи с нерегулярной структурой данных и продемонстрировать масштабируемость
на многопоточной системе с общей памятью (Shared Memory). Базовая математическая логика
(линейная комбинация столбцов) сохранена, изменена лишь схема распределения вычислений
между потоками.

## 2. Постановка задачи

- **Вход/Выход:** Аналогично SEQ (две матрицы CCS `A` и `B`, результат `C = A * B`).
- **Ограничения:** Совместимость размерностей проверяется в `ValidationImpl`.
- **Задача распараллеливания:** Безопасно распределить вычисления столбцов результирующей
  матрицы `C` между потоками, избежав состояния гонки (Data Race) при записи результатов.

## 3. Схема распараллеливания

В алгоритме умножения матриц в формате CCS (с подходом "столбец на столбец") вычисления
каждого столбца результирующей матрицы **абсолютно независимы** друг от друга.
Это делает внешний цикл по столбцам матрицы `B` идеальным кандидатом для распараллеливания.

### Выбранная область параллелизации и директивы

Параллелизуется внешний цикл по `cols_b_count`:

```cpp
#pragma omp parallel default(none) shared(mat_a, mat_b, thread_values, thread_row_indices, cols_b_count, rows_a_count)
```

- **`default(none)`**: Принуждает явно указывать статус каждой переменной (shared или private),
  что защищает от случайных гонок данных.
- **`shared`**: Матрицы `mat_a` и `mat_b` доступны только для чтения.
  Векторы `thread_values` и `thread_row_indices` заранее выделены по размеру `cols_b_count`,
  что позволяет потокам писать строго в "свои" ячейки без синхронизации.

### Изоляция памяти (Private Data)

Поскольку в последовательной версии использовался глобальный вектор `dense_col`,
в OMP-версии его пришлось локализовать. Каждый поток при входе в параллельную секцию
создает свои собственные буферы:

- `std::vector<std::complex<double>> accumulator`
- `std::vector<int> active_rows`
- `std::vector<int> row_marker`

Эти буферы переиспользуются потоком для всех столбцов, которые ему достанутся,
что минимизирует накладные расходы на аллокацию памяти.

### Планировщик (Schedule)

```cpp
#pragma omp for schedule(dynamic)
```

Использование `schedule(dynamic)` является критически важным решением.
В разреженных матрицах один столбец может быть пустым, а другой — плотно заполненным.
Статическое распределение (`static`) привело бы к сильному дисбалансу нагрузки
(Load Imbalance). Динамический планировщик выдает потокам новые столбцы
по мере их освобождения, обеспечивая равномерную загрузку ядер.

## 4. Детали реализации

```cpp
const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  auto &mat_res = GetOutput();

  const int rows_a_count = mat_a.count_rows;
  const int cols_b_count = mat_b.count_cols;

  std::vector<std::vector<std::complex<double>>> thread_values(static_cast<size_t>(cols_b_count));
  std::vector<std::vector<int>> thread_row_indices(static_cast<size_t>(cols_b_count));

#pragma omp parallel default(none) shared(mat_a, mat_b, thread_values, thread_row_indices, cols_b_count, rows_a_count)
  {
    std::vector<std::complex<double>> accumulator(static_cast<size_t>(rows_a_count), {0.0, 0.0});
    std::vector<int> active_rows;
    std::vector<int> row_marker(static_cast<size_t>(rows_a_count), -1);

#pragma omp for schedule(dynamic)
    for (int col_idx_b = 0; col_idx_b < cols_b_count; ++col_idx_b) {
      ProcessSingleColumn(col_idx_b, mat_a, mat_b, accumulator, active_rows, row_marker,
                          thread_values[static_cast<size_t>(col_idx_b)],
                          thread_row_indices[static_cast<size_t>(col_idx_b)]);
      active_rows.clear();
    }
  }

  mat_res.col_index[0] = 0;
  for (int col_idx = 0; col_idx < cols_b_count; ++col_idx) {
    const auto u_col_idx = static_cast<size_t>(col_idx);
    mat_res.values.insert(mat_res.values.end(), thread_values[u_col_idx].begin(), thread_values[u_col_idx].end());
    mat_res.row_index.insert(mat_res.row_index.end(), thread_row_indices[u_col_idx].begin(),
                             thread_row_indices[u_col_idx].end());
    mat_res.col_index[u_col_idx + 1] = static_cast<int>(mat_res.values.size());
  }

  return true;
```

## 5. Проверка корректности

Алгоритм является детерминированным на уровне содержимого каждого столбца.
Так как финальная склейка (конкатенация) происходит строго в последовательном порядке
от $0$ до $cols\_b\_count - 1$, результат многопоточной версии побитово совпадает
с последовательной (SEQ) для любых конфигураций потоков.
Код успешно прошел все функциональные тесты.

## 6. Экспериментальная среда

- **CPU:** Intel Core i5-12450H (8 ядер / 12 потоков)
- **RAM:** 16 GB
- **OS:** Ubuntu 24.04
- **Compiler:** GCC

## 7. Результаты производительности

Параметры теста: случайные разреженные матрицы 1000x1000 с плотностью заполнения 1%.

| Режим   | Потоки | Время Task_run (мс) | Ускорение | Эффективность |
| ------- | ------ | ------------------- | --------- | ------------- |
| SEQ     | 1      | 56.1                | 1.00x     | 100%          |
| OMP     | 2      | 31.8                | 1.76x     | 88.0%         |
| OMP     | 4      | 18.5                | 3.03x     | 75.7%         |

## 8. Выводы

1. **Хорошая масштабируемость:** На 2 потоках эффективность составила 88%,
   на 4 потоках — 75.7%.
2. **Динамическое планирование:** Использование `schedule(dynamic)` оправдало себя,
   позволив избежать простоя потоков на "пустых" столбцах.
3. **Узкие места:** Небольшое падение эффективности на 4 потоках связано с фазой
   постобработки. Несмотря на то, что вычисления распараллелены идеально,
   финальная конкатенация двумерных векторов `thread_values` в одномерный массив
   `mat_res.values` происходит последовательно в главном потоке и требует
   выделения памяти (реаллокации).

В целом, OpenMP-версия продемонстрировала высокую надежность и значительное сокращение
времени вычислений (ускорение в ~3 раза) при минимальном вмешательстве в архитектуру
последовательного кода.
