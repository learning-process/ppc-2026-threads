#include "ashihmin_d_mult_matr_crs/seq/include/ops_seq.hpp"

#include <cstddef>
#include <unordered_map>

namespace ashihmin_d_mult_matr_crs {

AshihminDMultMatrCrsSEQ::AshihminDMultMatrCrsSEQ(const InType &input_matrices) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_matrices;
}

bool AshihminDMultMatrCrsSEQ::ValidationImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;

  if (matrix_a.cols != matrix_b.rows) {
    return false;
  }
  if (matrix_a.row_ptr.size() != static_cast<std::size_t>(matrix_a.rows + 1)) {
    return false;
  }
  if (matrix_b.row_ptr.size() != static_cast<std::size_t>(matrix_b.rows + 1)) {
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

  GetOutput().rows = matrix_a.rows;
  GetOutput().cols = GetInput().second.cols;
  GetOutput().row_ptr.resize(matrix_a.rows + 1, 0);

  return true;
}

bool AshihminDMultMatrCrsSEQ::RunImpl() {
  const auto &matrix_a = GetInput().first;
  const auto &matrix_b = GetInput().second;
  auto &matrix_c = GetOutput();

  matrix_c.values.clear();
  matrix_c.col_index.clear();

  for (std::size_t row_index = 0; row_index < static_cast<std::size_t>(matrix_a.rows); ++row_index) {
    std::unordered_map<int, double> accumulator;

    for (std::size_t index_a = matrix_a.row_ptr[row_index];
         index_a < static_cast<std::size_t>(matrix_a.row_ptr[row_index + 1]); ++index_a) {
      int col_a = matrix_a.col_index[index_a];
      double value_a = matrix_a.values[index_a];

      for (std::size_t index_b = matrix_b.row_ptr[col_a];
           index_b < static_cast<std::size_t>(matrix_b.row_ptr[col_a + 1]); ++index_b) {
        int col_b = matrix_b.col_index[index_b];
        double value_b = matrix_b.values[index_b];
        accumulator[col_b] += value_a * value_b;
      }
    }

    for (const auto &[column_index, value] : accumulator) {
      if (value != 0.0) {
        matrix_c.values.push_back(value);
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
