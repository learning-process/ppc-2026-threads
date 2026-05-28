# Маркировка связанных компонент бинарного изображения — OMP

- Student: Кондрашова Виктория Андреевна
- Technology: OMP
- Variant: 29

## 1. Контекст

В данной реализации OpenMP используется для параллельного сканирования горизонтальных полос изображения.
В отличие от стандартного использования `omp for`, здесь применено ручное распределение нагрузки внутри параллельной области,
что дает полный контроль над декомпозицией данных.

## 2. Постановка задачи

Краткая постановка: маркировка связных областей. SEQ-версия ([seq/report.md](../seq/report.md)) — baseline.

## 3. Базовый алгоритм

Последовательная версия использует BFS. В OMP-версии изображение делится на $N$ полос, каждая из которых сканируется
своим потоком. Слияние результатов полос выполняется через DSU.

## 4. Схема распараллеливания

- **Область:** Параллельный блок `#pragma omp parallel`.
- **Директива:** `#pragma omp parallel num_threads(num_threads) default(none) ...`
- **Переменные:**
  - `shared`: `local_labels` (массив меток), `width`, `height`, `image`, `num_threads`.
  - `firstprivate`: `max_labels_per_thread` (смещение для уникальности меток в потоках).
- **Синхронизация:**
  - **Ручное распределение:** Потоки не используют автоматический `schedule`. Вместо этого каждый поток через
    `omp_get_thread_num()` вычисляет свой диапазон строк (`row_start`, `row_end`).
  - **Барьер:** Неявный барьер в конце области `parallel` гарантирует завершение всех сканирований перед этапом слияния.

## 5. Детали реализации

- **Файлы:** `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`
- **Особенности:** Для корректной работы в среде MSVC члены класса копируются в локальные переменные перед входом в
  параллельную область.
- **Гонки:** Отсутствуют, так как каждый поток пишет в свои строки массива `local_labels`.

```cpp
#pragma omp parallel num_threads(num_threads) default(none) shared(local_labels, width, height, image, num_threads) \
    firstprivate(max_labels_per_thread)
{
    const int tid = omp_get_thread_num();
    const int row_start = (tid * height) / num_threads;
    const int row_end = ((tid + 1) * height) / num_threads;
    const int label_offset = tid * max_labels_per_thread;
    ScanStripe(row_start, row_end, width, label_offset, image, local_labels);
}
```

## 6. Проверка корректности

Сравнение с SEQ на 13 тестах. Особое внимание тесту `boundary_bridge`, проверяющему корректность `MergeBoundaries`.

## 7. Экспериментальная среда

- **CPU:** AMD Ryzen 5 5600H (6 ядер).
- **RAM:** 16 GB.
- **OS:** Windows 11.
- **Compiler:** MSVC 19.37.
- **Build Type:** Release.
- **Команда:** `$env:PPC_NUM_THREADS=4; .\build\bin\ppc_perf_tests.exe --gtest_filter="*OMP*"`

## 8. Результаты

Замеры на `Chessboard` (512x512):

| Threads | Time (sec) | Speedup | Efficiency |
| :--- | :--- | :--- | :--- |
| 1 | 0.004414 | 1.00x | 100% |
| 2 | 0.003502 | 1.26x | 63.0% |
| 4 | 0.003079 | 1.43x | 35.8% |

## 9. Выводы

Ручное распределение нагрузки внутри OMP parallel секции обеспечило предсказуемое поведение алгоритма.
Основное ограничение масштабируемости — последовательный этап слияния строк через DSU.
