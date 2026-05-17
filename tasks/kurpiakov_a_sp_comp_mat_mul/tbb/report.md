**TBB report**  
# Умножение разреженных комплексных матриц в CSR — TBB  
- Student: Курпяков Алексей Георгиевич  
- Technology: TBB  
- Variant: 6  

## 1. Контекст  
TBB‑версия распараллеливает строки через `parallel_for`, используя потоковые буферы.

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
`parallel_for` по диапазону строк, буферы через `tbb::enumerable_thread_specific`.  
Grain size и partitioner не заданы явно.

Файл: tasks/kurpiakov_a_sp_comp_mat_mul/tbb/src/ops_tbb.cpp.  
Кодовый фрагмент:  
```cpp
tbb::enumerable_thread_specific<ThreadLocalRowBuffers> tls_buffers([&] { return ThreadLocalRowBuffers(cols); });

tbb::parallel_for(tbb::blocked_range<int>(0, rows), [&](const tbb::blocked_range<int> &range) {
  auto &buffers = tls_buffers.local();
  for (int i = range.begin(); i < range.end(); ++i) {
    ComputeRow(a, b, i, buffers, row_values[i], row_cols[i]);
  }
});
```

## 5. Детали реализации  
Сбор итоговых массивов выполняется параллельно в `BuildResultFromRows`.

## 6. Проверка корректности  
Сравнение с эталоном `Multiply` в perf‑тесте: tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp.

## 7. Экспериментальная среда  
Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`.

## 8. Результаты  
**pipeline**  
| Потоки | Время | Speedup | Efficiency |
|---:|---:|---:|---:|
| 4 | 7438 | 1.092 | 0.273 |

**task_run**  
| Потоки | Время | Speedup | Efficiency |
|---:|---:|---:|---:|
| 4 | 6843 | 1.243 | 0.311 |

## 9. Выводы  
TBB показал лучшее время среди измеренных, особенно в task_run.
