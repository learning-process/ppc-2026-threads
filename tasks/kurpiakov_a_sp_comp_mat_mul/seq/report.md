# Умножение разреженных комплексных матриц в формате CSR — SEQ

- Student: Курпяков Алексей Георгиевич  
- Technology: SEQ  
- Variant: 6  

## 1. Контекст

SEQ‑версия задает baseline корректности и времени для всех параллельных
реализаций. Реализация полностью последовательная, без параллельных участков.

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

## 3. Базовый алгоритм

Для каждой строки $i$ матрицы $A$ накапливаются вклады в строку результата:

1. Создаются буферы `row_acc` и `row_used` длиной `B.cols`.  
2. Для каждого ненулевого элемента $A(i, k)$ обходится строка $k$ матрицы $B$;
  сумма по колонкам накапливается в `row_acc`.  
3. Список использованных колонок сортируется и записывается в CSR‑структуру результата.

**Сложность:**  
$O\!\left(\sum_{i}\sum_{(i,k)\in A} nnz(B_k)\right) + O\!\left(\sum_i u_i \log u_i\right)$.  
**Память:** $O(B.cols)$ на буферы строки.

## 4. Детали реализации

Файлы: tasks/kurpiakov_a_sp_comp_mat_mul/seq/src/ops_seq.cpp,
tasks/kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp.

Кодовый фрагмент (валидация и запуск):

```cpp
bool KurpiskovACRSMatMulSEQ::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!ValidateCSR(a) || !ValidateCSR(b)) {
    return false;
  }

  return a.cols == b.rows;
}

bool KurpiskovACRSMatMulSEQ::RunImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = a.Multiply(b);
  return true;
}
```

## 5. Проверка корректности

Функциональные тесты SEQ: tasks/kurpiakov_a_sp_comp_mat_mul/tests/functional/main.cpp.  
Perf‑тесты сравнивают результат с эталоном `Multiply`: tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp.

## 6. Экспериментальная среда и методика

CPU: AMD Ryzen 7 7840HS, RAM: 4 GB, OS: Ubuntu 22.04, Compiler: GCC 13, build type: Release.  
`PPC_NUM_THREADS = 4` (для SEQ фактически используется 1 поток).  

**Определения метрик:**

- `time` — значение из perf‑логов (условные единицы, без пересчета).  
- $S = T_{seq}/T_{x}$ — ускорение для того же режима.  
- `workers` $W$ — число рабочих единиц (для SEQ $W = 1$).  
- $E = S / W$ — эффективность.

## 7. Результаты

| Mode | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| pipeline | 1 | 8125 | 1.000 | 1.000 |
| task_run | 1 | 8509 | 1.000 | 1.000 |

## 8. Выводы

SEQ‑версия — эталон без параллелизма и база для расчета ускорений в остальных отчетах.
