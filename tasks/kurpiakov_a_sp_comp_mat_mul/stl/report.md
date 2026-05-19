# Умножение разреженных комплексных матриц в CSR — STL

- Student: Курпяков Алексей Георгиевич  
- Technology: STL (`std::thread`)  
- Variant: 6  

## 1. Контекст

STL‑версия реализует ручное управление потоками и динамическое распределение строк.

## 2. Постановка задачи

См. SEQ: вход — две CSR‑матрицы с комплексными элементами, выход — их произведение в CSR.

## 3. Базовый алгоритм

См. SEQ: построчное накопление вкладов и сбор результата.

## 4. Схема распараллеливания

Динамическая выдача строк через `std::atomic<int> next_row`, локальные буферы
внутри потоков. `join()` вызывается после запуска всех потоков, что обеспечивает
реальный параллелизм.

Кодовый фрагмент (файл tasks/kurpiakov_a_sp_comp_mat_mul/stl/src/ops_stl.cpp):

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

## 7. Экспериментальная среда и методика

Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`.

**Определения метрик:**

- `time` — значение из perf‑логов (условные единицы).  
- $S = T_{seq}/T_{x}$ — ускорение для того же режима.  
- `workers` $W$ — число потоков STL.  
- $E = S / W$ — эффективность.

## 8. Результаты

| Mode | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: |
| pipeline | 4 | 9701 | 0.838 | 0.210 |
| task_run | 4 | 8545 | 0.996 | 0.249 |

## 9. Выводы

STL‑версия близка к SEQ; возможные причины — накладные расходы на создание
потоков и динамическую раздачу строк (без профилирования это гипотеза).
