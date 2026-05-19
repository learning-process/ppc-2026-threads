# Умножение разреженных комплексных матриц в CSR — ALL

- Student: Курпяков Алексей Георгиевич  
- Technology: ALL (MPI + std::thread)  
- Variant: 6  

## 1. Контекст

Гибридная версия: распределение строк между MPI‑процессами и многопоточная обработка внутри rank‑а.

## 2. Постановка задачи

См. SEQ: вход — две CSR‑матрицы с комплексными элементами, выход — их произведение в CSR.

## 3. Базовый алгоритм

См. SEQ: построчное накопление вкладов и сбор результата.

## 4. Межпроцессная схема

Ранги получают диапазоны строк. Для сборки общего CSR используются коллективные операции:

- `MPI_Allreduce` суммирует $nnz$ по строкам и синхронизирует процессы.  
- `MPI_Gatherv` собирает значения и индексы строк на rank 0.  
- `MPI_Bcast` распространяет собранный результат на все rank‑и.

Кодовый фрагмент (файл tasks/kurpiakov_a_sp_comp_mat_mul/all/src/ops_all.cpp):

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

## 7. Экспериментальная среда и методика

Ubuntu 22.04, GCC 13, Release, `PPC_NUM_THREADS = 4`, `PPC_NUM_PROC = 4`.

**Определения метрик:**

- `time` — значение из perf‑логов (условные единицы).  
- $S = T_{seq}/T_{x}$ — ускорение для того же режима.  
- `workers` $W$ — общее число рабочих единиц, $W = ranks \times threads$.  
- $E = S / W$ — эффективность.

## 8. Результаты

| Mode | Ranks | Threads/Rank | Workers | Time | Speedup | Efficiency |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| pipeline | 4 | 4 | 16 | 13759 | 0.590 | 0.037 |
| task_run | 4 | 4 | 16 | 12695 | 0.670 | 0.042 |

## 9. Выводы

ALL уступает SEQ на текущем наборе входов; наиболее вероятная причина —
стоимость коллективных коммуникаций и сборки глобального результата (без
профилирования это гипотеза).
