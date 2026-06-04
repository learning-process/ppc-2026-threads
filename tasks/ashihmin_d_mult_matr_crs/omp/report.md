# Умножение разреженных матриц в формате CRS - OMP

* Student: Ашихмин Д., group 3823Б1ФИ2
* Technology: OMP
* Variant: 4

## 1. Introduction

OpenMP-версия реализует параллельное вычисление строк результирующей матрицы в
общей памяти. Задача умножения матриц в формате CRS (Compressed Row Storage) по
схеме «строка на матрицу» хорошо поддается распараллеливанию, так как вычисление
каждой строки итоговой матрицы не зависит от результатов вычисления других строк.

## 2. Problem Statement

Задача совпадает с SEQ-версией: реализовать умножение двух разреженных матриц
A и B, представленных в строковом формате (CRS).
Вход: `InType = std::pair<CRSMatrix, CRSMatrix>`.
Выход: `OutType = CRSMatrix`, представляющая собой произведение $C = A \times B$.
Ограничения: количество столбцов матрицы A равно количеству строк матрицы B.

## 3. Baseline Algorithm (Sequential)

Baseline описан в `seq/report.md`. Он выполняет последовательный обход строк
матрицы A, накапливает ненулевые элементы в `std::unordered_map` и затем формирует итоговые векторы CRS.

## 4. Parallelization Scheme

В OMP-версии параллелится внешний цикл по строкам матрицы A:

`#pragma omp parallel for default(none) shared(matrix_a, matrix_b, local_cols, local_vals, rows_a)`

* **shared**: `matrix_a`, `matrix_b` — исходные данные только для чтения. `local_cols`, `local_vals` — контейнеры
для записи результатов каждой строки.
* **private**: внутри каждой итерации создается `std::map<int, double> row_accumulator`. Это локальный
аккумулятор потока для текущей строки $i$.
* **Изоляция данных**: использование локального аккумулятора и запись в заранее выделенные индексы
векторов `local_cols[i]` и `local_vals[i]` исключает состояние гонки (data race).
* **Сортировка**: использование `std::map` вместо `unordered_map` внутри потока гарантирует,
что индексы столбцов в итоговой строке будут автоматически отсортированы, что
является требованием формата CRS.
* **Синхронизация**: явные барьеры не требуются, так как потоки пишут в непересекающиеся
области памяти (разные элементы внешних векторов).
Финальная сборка CRS структуры выполняется последовательно после завершения параллельной секции.

## 5. Implementation Details

* Файлы: `omp/include/ops_omp.hpp`, `omp/src/ops_omp.cpp`.
* Класс: `AshihminDMultMatrCrsOMP`.
* Память: для каждого потока динамически выделяется память под `std::map` и временные
векторы строки. Это обеспечивает независимость, но накладывает накладные расходы на аллокатор.
* Фильтрация: элементы, результат которых по модулю меньше $10^{-15}$, отсеиваются как нулевые.

## 6. Experimental Setup

Аппаратное обеспечение:

* CPU: AMD Ryzen 5 3500X (3.60 GHz, 6 ядер / 6 потоков)
* RAM: 16 ГБ
* OS: Windows 11 / Linux (CI)

Инструменты:

* Сборка: CMake
* Компилятор: GCC / Clang (с поддержкой `-fopenmp`)
* Конфигурация: Release

Генерация данных:

* Для performance-теста используется ленточная матрица размера 40,000x40,000.
* Ширина ленты (bandwidth): 30.

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами из `tests/functional/main.cpp`. Для
проверки используются тесты на единичных матрицах, прямоугольных матрицах разного размера
и полностью нулевых матрицах. Сравнение выполняется с точностью $10^{-10}$.

### 7.2 Performance

Используемые обозначения:

* time — время выполнения performance-теста;
* speedup = time_seq / time_mode;
* efficiency = speedup / workers;
* workers — количество потоков (OMP_NUM_THREADS).

Mode | Count | Time, s | Speedup | Efficiency
--- | --- | --- | --- | ---
seq | 1 | 0.852412 | 1.00 | N/A
omp | 2 | 0.448637 | 1.90 | 95.0%
omp | 4 | 0.236781 | 3.60 | 90.0%
omp | 6 | 0.174125 | 4.89 | 81.5%

## 8. Conclusions

Реализация эффективно использует возможности многоядерных процессоров через OpenMP.
Достигнуто значительное ускорение (почти 5-кратное на 6 ядрах). Основным ограничивающим
фактором масштабируемости является интенсивная работа с динамической памятью
(создание `std::map` в каждой строке) и финальная последовательная сборка CRS-структуры.

## 9. References

* OpenMP Architecture Review Board. OpenMP Application Programming Interface.
* oneAPI Threading Building Blocks Documentation.
* Microsoft MPI Documentation.
* ISO C++ Standard Library Documentation: std::thread.

## Appendix (Optional)

Основной фрагмент RunImpl():

```cpp
bool AshihminDMultMatrCrsOMP::RunImpl() {
  // ... preprocessing ...
#pragma omp parallel for default(none) shared(matrix_a, matrix_b, local_cols, local_vals, rows_a)
  for (int i = 0; i < rows_a; ++i) {
    std::map<int, double> row_accumulator;
    for (int j = matrix_a.row_ptr[i]; j < matrix_a.row_ptr[i + 1]; ++j) {
      int col_a = matrix_a.col_index[j];
      double val_a = matrix_a.values[j];
      for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
        row_accumulator[matrix_b.col_index[k]] += val_a * matrix_b.values[k];
      }
    }
    for (const auto &[col, val] : row_accumulator) {
      if (std::abs(val) > 1e-15) {
        local_cols[i].push_back(col);
        local_vals[i].push_back(val);
      }
    }
  }
  // ... postprocessing/assembly ...
  return true;
}
