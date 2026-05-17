**STL report**  
# Умножение разреженных комплексных матриц в CSR — STL  
- Student: Курпяков Алексей Георгиевич  
- Technology: STL (`std::thread`)  
- Variant: 6  

## 1. Контекст  
STL‑версия реализует ручное управление потоками и динамическое распределение строк.

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
Динамическая выдача строк через `std::atomic<int> next_row`, локальные буферы внутри потоков, `join()` после запуска всех потоков.

Файл: tasks/kurpiakov_a_sp_comp_mat_mul/stl/src/ops_stl.cpp.  
Кодовый фрагмент:  
```cpp
std::atomic<int> next_row(0);
std::vector<std::thread> workers;

for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
  workers.emplace_back([&]() {
    ThreadLocalRowState state(b.cols);

    while (true) {
      const int row_idx = next_row.fetch_add(1, std::memory_order_relaxed);
      if (row_idx >= rows) {
        break;
      }

      MultiplySingleRow(a, b, row_idx, state, row_values[row_idx], row_cols[row_idx]);
    }
  });
}

for (auto &worker : workers) {
  worker.join();
}
```

## 5. Детали реализации  
Число потоков берется из `ppc::util::GetNumThreads()` и ограничивается числом строк.

## 6. Проверка корректности  
Сравнение с эталоном `Multiply` в perf‑тесте: tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp.

## 7. Экспериментальная среда  
Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`.

## 8. Результаты  
**pipeline**  
| Потоки | Время | Speedup | Efficiency |
|---:|---:|---:|---:|
| 4 | 9701 | 0.838 | 0.210 |

**task_run**  
| Потоки | Время | Speedup | Efficiency |
|---:|---:|---:|---:|
| 4 | 8545 | 0.996 | 0.249 |

## 9. Выводы  
STL‑версия близка к SEQ, но ускорение ограничено накладными расходами на управление потоками.
