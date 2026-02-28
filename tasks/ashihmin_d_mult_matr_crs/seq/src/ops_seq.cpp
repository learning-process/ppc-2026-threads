#include "ashihmin_d_mult_matr_crs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace ashihmin_d_mult_matr_crs {

AshihminDMultMatrCrsSEQ::AshihminDMultMatrCrsSEQ(const InType &input_data) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_data;
}

bool AshihminDMultMatrCrsSEQ::ValidationImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;

  if (matrix_a.cols != matrix_b.rows) {
    return false;
  }

  if (matrix_a.row_ptr.size() != static_cast<size_t>(matrix_a.rows + 1)) {
    return false;
  }

  if (matrix_b.row_ptr.size() != static_cast<size_t>(matrix_b.rows + 1)) {
    return false;
  }

  if (matrix_a.values.size() != matrix_a.col_index.size()) {
    return false;
  }

  if (matrix_b.values.size() != matrix_b.col_index.size()) {
    return false;
  }

  return true;
}

bool AshihminDMultMatrCrsSEQ::PreProcessingImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;

  auto &matrix_c = GetOutput();

  matrix_c.rows = matrix_a.rows;
  matrix_c.cols = matrix_b.cols;
  matrix_c.row_ptr.resize(matrix_a.rows + 1, 0);

  return true;
}

bool AshihminDMultMatrCrsSEQ::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();

  matrix_c.values.clear();
  matrix_c.col_index.clear();

  for (int row_index = 0; row_index < matrix_a.rows; ++row_index) {
    std::unordered_map<int, double> accumulator;

    int row_begin = matrix_a.row_ptr[row_index];
    int row_end = matrix_a.row_ptr[row_index + 1];

    for (int element_index = row_begin; element_index < row_end; ++element_index) {
      int column_a = matrix_a.col_index[element_index];
      double value_a = matrix_a.values[element_index];

      int b_row_begin = matrix_b.row_ptr[column_a];
      int b_row_end = matrix_b.row_ptr[column_a + 1];

      for (int b_element_index = b_row_begin; b_element_index < b_row_end; ++b_element_index) {
        int column_b = matrix_b.col_index[b_element_index];
        double value_b = matrix_b.values[b_element_index];

        accumulator[column_b] += value_a * value_b;
      }
    }

    std::vector<std::pair<int, double>> sorted_row(accumulator.begin(), accumulator.end());
    std::sort(sorted_row.begin(), sorted_row.end(),
              [](const auto &left, const auto &right) { return left.first < right.first; });

    for (const auto &[column_index, result_value] : sorted_row) {
      if (result_value != 0.0) {
        matrix_c.values.push_back(result_value);
        matrix_c.col_index.push_back(column_index);
      }
    }

    matrix_c.row_ptr[row_index + 1] = static_cast<int>(matrix_c.values.size());
  }

  return true;
}

bool AshihminDMultMatrCrsSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace ashihmin_d_mult_matr_crs
