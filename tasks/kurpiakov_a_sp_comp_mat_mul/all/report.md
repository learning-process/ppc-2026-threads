**ALL report**  
# Умножение разреженных комплексных матриц в CSR — ALL  
- Student: Курпяков Алексей Георгиевич  
- Technology: ALL (MPI + std::thread)  
- Variant: 6  

## 1. Контекст  
Гибридная версия: распределение строк между MPI‑процессами и многопоточная обработка внутри rank‑а.

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

## 4. Межпроцессная схема  
Ранги получают диапазоны строк, далее используют коллективные операции для сборки общего CSR.

Файл: tasks/kurpiakov_a_sp_comp_mat_mul/all/src/ops_all.cpp.  
Кодовый фрагмент:  
```cpp
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &world_size);

const auto [row_begin, row_end] = GetRowRange(rows, rank, world_size);

MPI_Allreduce(local_row_nnz.data(), global_row_nnz.data(), rows, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

MPI_Gatherv(local_re.data(), local_nnz, MPI_DOUBLE, global_re.data(), recv_counts.data(), recv_displs.data(),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

MPI_Bcast(global_re.data(), total_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
```

## 5. Внутрипроцессная схема  
Внутри rank‑а используется `std::thread` с динамической выдачей строк через `std::atomic<int> next_row`.

## 6. Проверка корректности  
Сравнение с эталоном `Multiply` в perf‑тесте: tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp.

## 7. Экспериментальная среда  
Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`, `PPC_NUM_PROC = 4`.

## 8. Результаты  
**pipeline**  
| Ранги | Потоки/ранг | Время | Speedup | Efficiency |
|---:|---:|---:|---:|---:|
| <PPC_NUM_PROC> | 4 | 13759 | 0.590 | 0.193 |

**task_run**  
| Ранги | Потоки/ранг | Время | Speedup | Efficiency |
|---:|---:|---:|---:|---:|
| <PPC_NUM_PROC> | 4 | 12695 | 0.670 | 0.193 |

## 9. Выводы  
ALL проигрывает SEQ из‑за стоимости коммуникаций и сборки полного результата у всех rank‑ов.

