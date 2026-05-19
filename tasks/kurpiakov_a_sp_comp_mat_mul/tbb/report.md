# Умножение разреженных комплексных матриц в CSR — TBB

- Student: Курпяков Алексей Георгиевич  
- Technology: TBB  
- Variant: 6  

## 1. Контекст

TBB‑версия распараллеливает строки через `parallel_for`, используя потоковые буферы.

## 2. Постановка задачи

См. SEQ: вход — две CSR‑матрицы с комплексными элементами, выход — их произведение в CSR.

## 3. Базовый алгоритм

См. SEQ: построчное накопление вкладов и сбор результата.

## 4. Схема распараллеливания

Диапазон строк $[0, rows)$ делится на поддиапазоны `tbb::blocked_range<int>`.
Каждый поддиапазон обрабатывается независимо, поэтому гонки отсутствуют.

**Grain size:** не задан явно, используется значение по умолчанию из TBB runtime.  
**Partitioner:** не задан, используется стандартный auto‑partitioner.  
**Контроль конкуренции:** `global_control` не используется; число worker‑потоков
определяется runtime (без явного ограничения в коде).

Кодовый фрагмент (файл tasks/kurpiakov_a_sp_comp_mat_mul/tbb/src/ops_tbb.cpp):

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

## 7. Экспериментальная среда и методика

Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`.

**Определения метрик:**

- `time` — значение из perf‑логов (условные единицы).  
- $S = T_{seq}/T_{x}$ — ускорение для того же режима.  
- `workers` $W$ — число потоков TBB.  
- $E = S / W$ — эффективность.

## 8. Результаты

| Mode | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| pipeline | 4 | 7438 | 1.092 | 0.273 |
| task_run | 4 | 6843 | 1.243 | 0.311 |

## 9. Выводы

TBB показывает наилучшее время среди измеренных; вывод основан на приведенных таблицах без дополнительного профилирования.
