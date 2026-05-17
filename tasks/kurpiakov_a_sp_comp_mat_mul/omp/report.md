
**OMP report**  
# Умножение разреженных комплексных матриц в CSR — OMP  
- Student: Курпяков Алексей Георгиевич  
- Technology: OMP  
- Variant: 6  

## 1. Контекст  
OMP‑версия распараллеливает независимые строки результата, сохраняя алгоритм SEQ.

## 2. Постановка задачи  
**Вход:** пара матриц $A$ и $B$ в формате CSR с элементами типа `Complex<double>`.  
**Выход:** матрица $C = A \times B$ в формате CSR.

**Ограничения и корректность входа:**
- $rows > 0$, $cols > 0$ для каждой матрицы.  
- `row_ptr` имеет длину $rows + 1$, `row_ptr[0] = 0`.  
- `row_ptr[rows]` равен `values.size()`.  
- `col_indices.size()` равен `values.size()`.  
- Каждый `col_indices[j]` лежит в диапазоне $[0, cols)$.  
- Размерности должны удовлетворять $A.cols = B.rows$.  
- Неявно предполагается монотонность `row_ptr` (иначе обход CSR некорректен).

## 3. Базовый алгоритм  
Для каждой строки $i$ матрицы $A$ алгоритм накапливает вклад в строку результата:

1. Создаются буферы `row_acc` (аккумулятор значений по всем колонкам $B$) и `row_used` (отметки использованных колонок).  
2. Для каждого ненулевого элемента $A(i, k)$ обходится строка $k$ матрицы $B$; вклад $A(i, k) \cdot B(k, j)$ накапливается в `row_acc[j]`.  
3. Список использованных колонок сортируется и записывается в CSR‑структуру результата.

**Сложность:**  
$O\!\left(\sum_{i}\sum_{(i,k)\in A} nnz(B_k)\right) + O\!\left(\sum_i u_i \log u_i\right)$, где $u_i$ — число затронутых колонок в строке $i$.  
**Память:** $O(B.cols)$ на `row_acc` и `row_used`.

## 4. Схема распараллеливания  
Параллелизация по строкам; каждая строка заполняется одним потоком.  
Используется `#pragma omp parallel` и `#pragma omp for schedule(dynamic)`.

Файл: tasks/kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp.  
Кодовый фрагмент:  
```cpp
#pragma omp parallel default(none) shared(self, other, row_values, row_col_indices, nrows, ncols)
{
  std::vector<T> acc_re(static_cast<std::size_t>(ncols));
  std::vector<T> acc_im(static_cast<std::size_t>(ncols));
  std::vector<bool> local_used(static_cast<std::size_t>(ncols), false);

#pragma omp for schedule(dynamic)
  for (int i = 0; i < nrows; ++i) {
    ProcessRow(i, self, other, acc_re, acc_im, local_used, row_values, row_col_indices);
  }
}
```

## 5. Детали реализации  
`RunImpl` вызывает `a.OMPMultiply(b)` в tasks/kurpiakov_a_sp_comp_mat_mul/omp/src/ops_omp.cpp.

## 6. Проверка корректности  
В perf‑тесте результат сравнивается с эталоном `Multiply`: tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp.

## 7. Экспериментальная среда  
Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`.

## 8. Результаты  
**pipeline**  
| Потоки | Время | Speedup | Efficiency |
|---:|---:|---:|---:|
| 4 | 9358 | 0.868 | 0.217 |

**task_run**  
| Потоки | Время | Speedup | Efficiency |
|---:|---:|---:|---:|
| 4 | 9036 | 0.942 | 0.236 |

## 9. Выводы  
На выбранном размере данных OMP проигрывает SEQ, что указывает на высокие накладные расходы и ограничение по памяти.