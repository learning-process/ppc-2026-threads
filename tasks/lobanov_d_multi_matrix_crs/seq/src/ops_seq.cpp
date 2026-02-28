#include "lobanov_d_multi_matrix_crs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"

namespace lobanov_d_multi_matrix_crs {

constexpr double kEpsilonThreshold = 1e-12;

LobanovMultyMatrixSEQ::LobanovMultyMatrixSEQ(const InType &input_matrices) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_matrices;
  CompressedRowMatrix empty_matrix;
  empty_matrix.row_count = 0;
  empty_matrix.column_count = 0;
  empty_matrix.non_zero_count = 0;
  GetOutput() = empty_matrix;
}

bool LobanovMultyMatrixSEQ::ValidationImpl() {
  const auto &[matrix_a, matrix_b] = GetInput();
  return (matrix_a.row_count > 0 && matrix_a.column_count > 0 && matrix_b.row_count > 0 && matrix_b.column_count > 0 &&
          matrix_a.column_count == matrix_b.row_count);
}

bool LobanovMultyMatrixSEQ::PreProcessingImpl() {
  return true;
}

namespace {

void MultiplyRowByMatrix(const std::vector<double> &a_row_values, const std::vector<int> &a_row_columns,
                         const CompressedRowMatrix &matrix_b, std::vector<double> &temp_row,
                         std::vector<int> &temp_col_markers, int row_index, std::vector<double> &result_values,
                         std::vector<int> &result_col_indices, int &result_row_start) {
  // Обнуляем временный массив для этой строки
  for (size_t i = 0; i < temp_row.size(); ++i) {
    if (temp_col_markers[i] == row_index) {
      temp_row[i] = 0.0;
      temp_col_markers[i] = -1;
    }
  }

  // Умножаем строку A на матрицу B
  for (size_t k = 0; k < a_row_columns.size(); ++k) {
    int b_row = a_row_columns[k];
    double a_value = a_row_values[k];

    for (int b_ptr = matrix_b.row_pointer_data[b_row]; b_ptr < matrix_b.row_pointer_data[b_row + 1]; ++b_ptr) {
      int b_col = matrix_b.column_index_data[b_ptr];
      double b_value = matrix_b.value_data[b_ptr];

      if (temp_col_markers[b_col] != row_index) {
        temp_col_markers[b_col] = row_index;
        temp_row[b_col] = a_value * b_value;
      } else {
        temp_row[b_col] += a_value * b_value;
      }
    }
  }

  // Собираем ненулевые элементы результата
  for (int col = 0; col < matrix_b.column_count; ++col) {
    if (temp_col_markers[col] == row_index && std::abs(temp_row[col]) > kEpsilonThreshold) {
      result_values.push_back(temp_row[col]);
      result_col_indices.push_back(col);
    }
  }

  result_row_start = static_cast<int>(result_values.size());
}

}  // namespace

void LobanovMultyMatrixSEQ::PerformMatrixMultiplication(const CompressedRowMatrix &first_matrix,
                                                        const CompressedRowMatrix &second_matrix,
                                                        CompressedRowMatrix &product_result) {
  product_result.row_count = first_matrix.row_count;
  product_result.column_count = second_matrix.column_count;

  product_result.row_pointer_data.clear();
  product_result.row_pointer_data.push_back(0);

  product_result.value_data.clear();
  product_result.column_index_data.clear();

  // Временные массивы для вычисления одной строки результата
  int cols = second_matrix.column_count;
  std::vector<double> temp_row(cols, 0.0);
  std::vector<int> temp_col_markers(cols, -1);

  // Обрабатываем каждую строку первой матрицы
  for (int i = 0; i < first_matrix.row_count; ++i) {
    int row_start = first_matrix.row_pointer_data[i];
    int row_end = first_matrix.row_pointer_data[i + 1];

    if (row_start == row_end) {
      // Пустая строка в A
      product_result.row_pointer_data.push_back(static_cast<int>(product_result.value_data.size()));
      continue;
    }

    // Извлекаем значения и индексы колонок для текущей строки A
    std::vector<double> a_row_values;
    std::vector<int> a_row_columns;

    for (int ptr = row_start; ptr < row_end; ++ptr) {
      a_row_values.push_back(first_matrix.value_data[ptr]);
      a_row_columns.push_back(first_matrix.column_index_data[ptr]);
    }

    // Умножаем строку на матрицу B
    int result_row_start = 0;
    MultiplyRowByMatrix(a_row_values, a_row_columns, second_matrix, temp_row, temp_col_markers, i,
                        product_result.value_data, product_result.column_index_data, result_row_start);

    product_result.row_pointer_data.push_back(static_cast<int>(product_result.value_data.size()));
  }

  product_result.non_zero_count = static_cast<int>(product_result.value_data.size());
}

bool LobanovMultyMatrixSEQ::RunImpl() {
  const auto &[matrix_a, matrix_b] = GetInput();

  try {
    CompressedRowMatrix result_matrix;
    result_matrix.row_count = 0;
    result_matrix.column_count = 0;
    result_matrix.non_zero_count = 0;

    PerformMatrixMultiplication(matrix_a, matrix_b, result_matrix);
    GetOutput() = result_matrix;
    return true;
  } catch (const std::exception &) {
    return false;
  }
}

bool LobanovMultyMatrixSEQ::PostProcessingImpl() {
  const auto &result_matrix = GetOutput();
  return result_matrix.row_count > 0 && result_matrix.column_count > 0 &&
         result_matrix.row_pointer_data.size() == static_cast<size_t>(result_matrix.row_count) + 1;
}

}  // namespace lobanov_d_multi_matrix_crs
