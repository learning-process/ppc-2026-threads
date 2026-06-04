# Умножение разреженных матриц в формате CRS (Compressed Row Storage)

* Student: Ашихмин Д., group 3823Б1ФИ2
* Technology: SEQ, OMP, TBB, STL, ALL
* Variant: 4

## 1. Introduction

Работа посвящена реализации высокопроизводительного алгоритма умножения разреженных матриц,
хранящихся в строковом формате CRS. Этот формат является стандартом для работы с большими
разреженными системами, так как позволяет хранить только ненулевые элементы. В работе реализованы
пять вариантов задачи: последовательный baseline и параллельные версии с использованием
OpenMP, TBB, STL и гибридной схемы ALL (MPI + Threads).

Ожидаемый результат работы — корректное произведение матриц и сравнение эффективности различных
технологий параллелизма при масштабировании на 6-ядерном процессоре.

## 2. Problem Statement

Входные данные задаются типом `InType = std::pair<CRSMatrix, CRSMatrix>`:

* matrix_a — левый операнд;
* matrix_b — правый операнд.

Выходные данные: `OutType = CRSMatrix`, результат операции $C = A \times B$.

Формат CRS включает три основных вектора: `values` (значения), `col_index` (номера столбцов)
и `row_ptr` (индексы начала строк). Ограничения: число столбцов матрицы A должно быть
равно числу строк матрицы B (`matrix_a.cols == matrix_b.rows`). Элементы имеют тип `double`.

## 3. Baseline Algorithm (Sequential)

Последовательная версия служит baseline для всех параллельных реализаций. Алгоритм реализует
схему «строка на матрицу». Для каждой строки $i$ матрицы A:

1. Используется временный аккумулятор `std::unordered_map<int, double>` для накопления значений.
2. Просматриваются ненулевые элементы $A_{ik}$.
3. Для каждого $A_{ik}$ выбирается соответствующая строка $k$ матрицы B.
4. Вычисляются произведения $A_{ik} \times B_{kj}$ и суммируются в аккумуляторе по индексу $j$.

В конце обработки строки элементы аккумулятора сортируются по индексу столбца и переносятся
в итоговую структуру. Время выполнения SEQ считается базовым (speedup = 1.0).

## 4. Parallelization Scheme

* **SEQ**: один поток, параллелизм не используется.
* **OMP**: параллельный цикл по строкам матрицы A с использованием `#pragma omp parallel for`.
  Каждый поток имеет локальный аккумулятор `std::map`.
* **TBB**: использование `tbb::parallel_for` для автоматического распределения строк между
  потоками планировщиком TBB (work-stealing).
* **STL**: ручное разбиение диапазона строк на блоки (chunks) и запуск их через `std::async`
  с политикой `launch::async`.
* **ALL**: гибридная схема. MPI делит строки между процессами. Внутри процесса используются
  STL-потоки для разделения задач, TBB для мелкозернистого параллелизма и OpenMP для обработки
  метаданных. Глобальная сборка разреженной матрицы выполняется через `MPI_Allgatherv`.

## 5. Implementation Details

Код задачи расположен в папке `tasks/ashihmin_d_mult_matr_crs`:

* `common/include/common.hpp` — структуры данных и типы;
* `seq/src/ops_seq.cpp` — последовательный baseline;
* `omp/src/ops_omp.cpp` — OpenMP реализация;
* `tbb/src/ops_tbb.cpp` — TBB реализация;
* `stl/src/ops_stl.cpp` — STL (std::async) реализация;
* `all/src/ops_all.cpp` — гибридная версия ALL;
* `tests/` — тесты производительности и корректности.

Во всех параллельных версиях вместо `unordered_map` применен `std::map`, что гарантирует
автоматическую сортировку столбцов согласно спецификации CRS.

## 6. Experimental Setup

Аппаратное обеспечение:

* CPU: AMD Ryzen 5 3500X (3.60 GHz, 6 ядер / 6 потоков)
* RAM: 16 ГБ
* OS: Windows 11 Pro / Linux (CI)

Инструменты:

* Сборка: CMake
* Компилятор: GCC / Clang / MSVC
* Конфигурация: Release

Окружение:

* `PPC_NUM_THREADS`: задает число потоков (1–6).
* `PPC_NUM_PROC`: задает число MPI-процессов для ALL.

Генерация данных:

* Используются ленточные матрицы размера 40,000x40,000 с шириной полосы 30.

## 7. Results and Discussion

### 7.1 Correctness

Корректность проверялась функциональными тестами. Вычисленное произведение сравнивается с
эталонным результатом (полученным последовательно) для единичных, прямоугольных и разреженных
ленточных матриц. Допустимая погрешность: $10^{-10}$.

### 7.2 Performance

Mode | Count | Time, s | Speedup | Efficiency
--- | --- | --- | --- | ---
seq | 1 | 0.852412 | 1.00 | N/A
omp | 6 | 0.174125 | 4.89 | 81.50%
tbb | 6 | 0.165214 | 5.16 | 86.00%
stl | 6 | 0.181245 | 4.70 | 78.33%
all | 2 x 3 | 0.192451 | 4.43 | 73.83%

## 8. Conclusions

Все параллельные реализации продемонстрировали существенное ускорение. Наилучший результат
показала технология TBB за счет балансировки нагрузки. Основное ограничение масштабируемости
связано с накладными расходами на выделение памяти для локальных аккумуляторов строк и
финальную последовательную сборку CRS-структуры в общей памяти.

## 9. References

* OpenMP Architecture Review Board. OpenMP Application Programming Interface.
* oneAPI Threading Building Blocks Documentation.
* Microsoft MPI Documentation.
* ISO C++ Standard Library Documentation: std::thread, std::async.

## Appendix (Optional)

Ниже приведены основные фрагменты RunImpl() для всех реализаций.

### SEQ RunImpl

```cpp
for (int row_index = 0; row_index < matrix_a.rows; ++row_index) {
  std::unordered_map<int, double> accumulator;
  for (int j = matrix_a.row_ptr[row_index]; j < matrix_a.row_ptr[row_index + 1]; ++j) {
    int col_a = matrix_a.col_index[j];
    double val_a = matrix_a.values[j];
    for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
      accumulator[matrix_b.col_index[k]] += val_a * matrix_b.values[k];
    }
  }
}
```

### OMP RunImpl

```cpp
#pragma omp parallel for default(none) shared(matrix_a, matrix_b, local_cols, local_vals, rows_a)
for (int i = 0; i < rows_a; ++i) {
  std::map<int, double> row_accumulator;
  for (int j = matrix_a.row_ptr[i]; j < matrix_a.row_ptr[i + 1]; ++j) {
    int col_a = matrix_a.col_index[j];
    for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
      row_accumulator[matrix_b.col_index[k]] += matrix_a.values[j] * matrix_b.values[k];
    }
  }
}
```

### TBB RunImpl

```cpp
tbb::parallel_for(0, rows_a, [&](int i) {
  std::map<int, double> row_accumulator;
  for (int j = matrix_a.row_ptr[i]; j < matrix_a.row_ptr[i + 1]; ++j) {
    int col_a = matrix_a.col_index[j];
    for (int k = matrix_b.row_ptr[col_a]; k < matrix_b.row_ptr[col_a + 1]; ++k) {
      row_accumulator[matrix_b.col_index[k]] += matrix_a.values[j] * matrix_b.values[k];
    }
  }
});
```

### STL RunImpl

```cpp
for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
  futures.push_back(std::async(std::launch::async, [=, &matrix_a, &matrix_b, &local_cols, &local_vals] {
    for (int i = start_row; i < end_row; ++i) {
      MultiplyRow(i, matrix_a, matrix_b, local_cols[i], local_vals[i]);
    }
  }));
}
```

### ALL RunImpl

```cpp
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
// Process local stripes of matrix A
for (int t = 0; t < thread_count; ++t) {
  threads.emplace_back(compute_rows, s, e); // Inside: tbb::parallel_for
}
for (auto &th : threads) th.join();
// Collect results using MPI_Allgatherv
MPI_Allgatherv(my_flat_cols.data(), ..., matrix_c.col_index.data(), ...);
