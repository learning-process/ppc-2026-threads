# Умножение разреженных матриц в формате CRS - SEQ

* Student: Ашихмин Д., group 3823Б1ФИ2
* Technology: SEQ
* Variant: 4

## 1. Introduction
Последовательная версия используется как базовая линия производительности. Она реализует классический алгоритм умножения разреженных матриц в формате CRS, который служит эталоном для проверки корректности результатов параллельных реализаций и расчета ускорения.

## 2. Problem Statement
Вход: `InType = std::pair<CRSMatrix, CRSMatrix>`, где каждая матрица представлена структурой с полями `rows`, `cols`, `row_ptr`, `col_index`, `values`.
Выход: `OutType = CRSMatrix`, результат умножения матриц A и B.
Ограничения: `matrix_a.cols == matrix_b.rows`, элементы типа `double`.

## 3. Baseline Algorithm (Sequential)
SEQ-версия использует алгоритм «строка на матрицу». Для каждой строки $i$ матрицы A создается временный аккумулятор (`std::unordered_map`). Алгоритм проходит по всем ненулевым элементам строки $i$ матрицы A. Для каждого элемента $A_{ik}$ находится соответствующая строка $k$ матрицы B, и её значения, умноженные на $A_{ik}$, добавляются в аккумулятор. 
После завершения обработки строки, данные из аккумулятора сортируются по индексам столбцов, фильтруются от нулевых значений и записываются в итоговую структуру CRS.

## 4. Parallelization Scheme
Параллелизм отсутствует. Все операции выполняются последовательно в одном потоке. `workers = 1`.

## 5. Implementation Details
* Файлы: `seq/include/ops_seq.hpp`, `seq/src/ops_seq.cpp`.
* Класс: `AshihminDMultMatrCrsSEQ`.
* `ValidationImpl()`: Проверяет совместимость матриц по размерности.
* `RunImpl()`: Реализует основной расчет с использованием `std::unordered_map` для накопления строк и `std::ranges::sort` для упорядочивания столбцов.

## 6. Experimental Setup
Аппаратное обеспечение:
* CPU: 12th Gen Intel(R) Core(TM) i5-12450H (8 ядер / 12 потоков)
* RAM: 16 ГБ
* OS: Windows 11 / Linux (CI)

Генерация данных:
* Для performance-теста используется ленточная матрица размера 40,000x40,000.
* Ширина ленты (bandwidth): 30.

## 7. Results and Discussion
### 7.1 Correctness
Корректность проверялась функциональными тестами. Вычисленные значения сравниваются с эталонными результатами для единичных и прямоугольных матриц с точностью $10^{-10}$.

### 7.2 Performance
Mode | Count | Time, s | Speedup | Efficiency
--- | --- | --- | --- | ---
seq | 1 | 0.852412 | 1.00 | N/A

## 8. Conclusions
Последовательная версия обеспечивает корректное перемножение разреженных матриц. Она является baseline для оценки эффективности многопоточных версий.

## 9. References
1. Microsoft MPI Documentation.
2. ISO C++ Standard Library Documentation: std::thread.
3. Sparse Matrix-Matrix Multiplication Algorithms.

## Appendix (Optional)
Основной фрагмент RunImpl():
```cpp
bool AshihminDMultMatrCrsSEQ::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();

  matrix_c.values.clear();
  matrix_c.col_index.clear();

  for (int row_index = 0; row_index < matrix_a.rows; ++row_index) {
    std::unordered_map<int, double> accumulator;

    auto row_begin = static_cast<std::size_t>(matrix_a.row_ptr[row_index]);
    auto row_end = static_cast<std::size_t>(matrix_a.row_ptr[row_index + 1]);

    for (std::size_t index_a = row_begin; index_a < row_end; ++index_a) {
      int col_a = matrix_a.col_index[index_a];
      double value_a = matrix_a.values[index_a];

      auto col_begin = static_cast<std::size_t>(matrix_b.row_ptr[col_a]);
      auto col_end = static_cast<std::size_t>(matrix_b.row_ptr[col_a + 1]);

      for (std::size_t index_b = col_begin; index_b < col_end; ++index_b) {
        int col_b = matrix_b.col_index[index_b];
        double value_b = matrix_b.values[index_b];

        accumulator[col_b] += value_a * value_b;
      }
    }

    std::vector<std::pair<int, double>> sorted_values(accumulator.begin(), accumulator.end());

    std::ranges::sort(sorted_values, {}, &std::pair<int, double>::first);

    for (const auto &entry : sorted_values) {
      if (entry.second != 0.0) {
        matrix_c.values.push_back(entry.second);
        matrix_c.col_index.push_back(entry.first);
      }
    }

    matrix_c.row_ptr[row_index + 1] = static_cast<int>(matrix_c.values.size());
  }

  return true;
}