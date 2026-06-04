
# Умножение разреженных комплексных матриц в CSR — OMP

- Student: Курпяков Алексей Георгиевич  
- Technology: OMP  
- Variant: 6  

## 1. Контекст

OMP‑версия распараллеливает независимые строки результата, сохраняя алгоритм SEQ.

## 2. Постановка задачи

См. SEQ: вход — две CSR‑матрицы с комплексными элементами, выход — их произведение в CSR.

## 3. Базовый алгоритм

См. SEQ: построчное накопление вкладов и сбор результата в CSR.

## 4. Схема распараллеливания

Параллелизация идет по строкам результата; каждая строка обрабатывается одним потоком.

**Атрибуты переменных и директивы:**

- `shared`: `self`, `other`, `row_values`, `row_col_indices`, `nrows`, `ncols`.  
- `private`: `acc_re`, `acc_im`, `local_used`, индекс `i` цикла.  
- `reduction`: не используется, так как каждая строка записывается в собственный буфер.  
- `schedule`: `dynamic`.

Кодовый фрагмент (файл tasks/kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp):

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

Результат сравнивается с эталоном `Multiply` в perf‑тесте: tasks/kurpiakov_a_sp_comp_mat_mul/tests/performance/main.cpp.

## 7. Экспериментальная среда и методика

Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`.

**Определения метрик:**

- `time` — значение из perf‑логов (условные единицы).  
- $S = T_{seq}/T_{x}$ — ускорение для того же режима.  
- `workers` $W$ — число потоков (OMP).  
- $E = S / W$ — эффективность.

## 8. Результаты

| Mode | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| pipeline | 4 | 9358 | 0.868 | 0.217 |
| task_run | 4 | 9036 | 0.942 | 0.236 |

## 9. Выводы

OMP‑результат ниже SEQ на выбранном размере; это согласуется с высокими
накладными расходами при построчной динамической балансировке без
подтверждающего профилирования.
