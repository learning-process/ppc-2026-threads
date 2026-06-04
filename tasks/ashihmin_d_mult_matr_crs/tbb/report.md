# Умножение разреженных матриц в формате CRS - TBB

* Student: Ашихмин Д., group 3823Б1ФИ2
* Technology: TBB
* Variant: 4

## 1. Introduction
TBB-версия использует задачно-ориентированный подход Intel oneAPI Threading Building Blocks. Вместо явного создания потоков, общий диапазон строк матрицы передаётся планировщику TBB, который разбивает его на подзадачи и распределяет работу между потоками, используя алгоритм кражи задач (work-stealing).

## 2. Problem Statement
Требуется реализовать умножение двух разреженных матриц A и B, представленных в строковом формате CRS (Compressed Row Storage). 
Вход: InType = std::pair<CRSMatrix, CRSMatrix>.
Выход: OutType = CRSMatrix, произведение C = A * B.
Ограничения: matrix_a.cols == matrix_b.rows, элементы типа double.

## 3. Baseline Algorithm (Sequential)
Baseline описан в seq/report.md. Один поток последовательно обходит строки матрицы A, накапливает ненулевые значения во временном словаре и формирует результирующие векторы values и col_index для каждой строки.

## 4. Parallelization Scheme
Используется функция tbb::parallel_for для распределения вычислений строк результирующей матрицы.

* **Диапазон:** tbb::parallel_for принимает диапазон индексов от 0 до количества строк матрицы A.
* **Планировщик:** Используется стандартный механизм TBB. Это эффективно для разреженных матриц, так как плотность строк может сильно различаться, а динамическая балансировка нагрузки (work-stealing) позволяет избежать простоя ядер.
* **Локальность данных:** Каждый рабочий поток внутри лямбда-выражения создает свой std::map<int, double> row_accumulator. Это обеспечивает потокобезопасность без использования блокировок.
* **Автоматическая сортировка:** Использование std::map гарантирует, что индексы столбцов внутри каждой строки будут отсортированы, что необходимо для корректности формата CRS.
* **Контроль конкуренции:** Количество используемых потоков ограничивается через настройки тестового фреймворка PPC (переменная PPC_NUM_THREADS).

## 5. Implementation Details
* Файлы: tbb/include/ops_tbb.hpp, tbb/src/ops_tbb.cpp.
* Класс: AshihminDMultMatrCrsTBB.
* Результаты каждой итерации записываются в заранее подготовленные локальные векторы local_cols[i] и local_vals[i].
* Гонок данных нет, так как каждый таск работает с уникальным индексом строки i и пишет в свою область памяти.
* Сборка финальной структуры (объединение локальных векторов в один плоский массив) выполняется последовательно после завершения параллельного цикла.

## 6. Experimental Setup
Аппаратное обеспечение:
* CPU: AMD Ryzen 5 3500X (3.60 GHz, 6 ядер)
* RAM: 16 ГБ
* OS: Windows 11 / Linux (CI)

Инструменты:
* Сборка: CMake
* Компилятор: GCC / Clang / MSVC
* Библиотека: Intel oneAPI TBB

## 7. Results and Discussion
7.1 Correctness
Корректность проверялась функциональными тестами. Вычисленные значения сравниваются с результатами последовательной версии для различных конфигураций: единичные, диагональные и прямоугольные матрицы.

7.2 Performance
Mode | Count | Time, s | Speedup | Efficiency
--- | --- | --- | --- | ---
seq | 1 | 0.852412 | 1.00 | N/A
tbb | 2 | 0.441664 | 1.93 | 96.5%
tbb | 4 | 0.226705 | 3.76 | 94.0%
tbb | 6 | 0.165214 | 5.16 | 86.0%

## 8. Conclusions
Реализация на TBB показала высокую эффективность и отличную масштабируемость. За счет механизма динамического распределения задач (work-stealing), версия TBB эффективно справляется с неравномерной плотностью строк в разреженных матрицах.

## 9. References
1. oneAPI Threading Building Blocks Documentation.
2. OpenMP Architecture Review Board. OpenMP API.
3. Microsoft MPI Documentation.
4. ISO C++ Standard Library Documentation: std::thread.

## Appendix (Optional)
Основной фрагмент RunImpl():
```cpp
tbb::parallel_for(0, rows_a, [&](int i) {
    std::map<int, double> row_accumulator;
    for (int j = matrix_a.row_ptr[i]; j < matrix_a.row_ptr[i + 1]; ++j) {
      int col_a = matrix_a.col_index[j];
      double val_a = matrix_a.values[j];
      for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
        row_accumulator[matrix_b.col_index[k]] += val_a * matrix_b.values[k];
      }
    }
    // Filter and save...
});