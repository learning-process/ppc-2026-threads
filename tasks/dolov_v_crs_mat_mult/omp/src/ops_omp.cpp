#include "dolov_v_crs_mat_mult/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <utility>
#include <vector>

#include "dolov_v_crs_mat_mult/common/include/common.hpp"

namespace dolov_v_crs_mat_mult {

DolovVCrsMatMultOmp::DolovVCrsMatMultOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DolovVCrsMatMultOmp::ValidationImpl() {
  const auto &input_data = GetInput();
  if (input_data.size() != 2) {
    return false;
  }
  const auto &matrix_a = input_data[0];
  const auto &matrix_b = input_data[1];
  return matrix_a.num_cols == matrix_b.num_rows && matrix_a.num_rows > 0 && matrix_b.num_cols > 0;
}

bool DolovVCrsMatMultOmp::PreProcessingImpl() {
  return true;
}

SparseMatrix DolovVCrsMatMultOmp::TransposeMatrix(const SparseMatrix &matrix) {
  SparseMatrix transposed;
  transposed.num_rows = matrix.num_cols;
  transposed.num_cols = matrix.num_rows;
  transposed.row_pointers.assign(transposed.num_rows + 1, 0);

  for (int col_idx : matrix.col_indices) {
    transposed.row_pointers[col_idx + 1]++;
  }
  for (int i = 0; i < transposed.num_rows; ++i) {
    transposed.row_pointers[i + 1] += transposed.row_pointers[i];
  }

  transposed.values.resize(matrix.values.size());
  transposed.col_indices.resize(matrix.col_indices.size());

  std::vector<int> current_pos = transposed.row_pointers;
  for (int i = 0; i < matrix.num_rows; ++i) {
    for (int j = matrix.row_pointers[i]; j < matrix.row_pointers[i + 1]; ++j) {
      int col = matrix.col_indices[j];
      int dest_idx = current_pos[col]++;
      transposed.values[dest_idx] = matrix.values[j];
      transposed.col_indices[dest_idx] = i;
    }
  }
  return transposed;
}

double DolovVCrsMatMultOmp::DotProduct(const SparseMatrix &matrix_a, int row_a, const SparseMatrix &matrix_b_t,
                                       int row_b) {
  double sum = 0.0;
  int ptr_a = matrix_a.row_pointers[row_a];
  int ptr_b = matrix_b_t.row_pointers[row_b];
  const int end_a = matrix_a.row_pointers[row_a + 1];
  const int end_b = matrix_b_t.row_pointers[row_b + 1];

  while (ptr_a < end_a && ptr_b < end_b) {
    if (matrix_a.col_indices[ptr_a] == matrix_b_t.col_indices[ptr_b]) {
      sum += matrix_a.values[ptr_a] * matrix_b_t.values[ptr_b];
      ptr_a++;
      ptr_b++;
    } else if (matrix_a.col_indices[ptr_a] < matrix_b_t.col_indices[ptr_b]) {
      ptr_a++;
    } else {
      ptr_b++;
    }
  }
  return sum;
}

bool DolovVCrsMatMultOmp::RunImpl() {
  const auto &matrix_a = GetInput()[0];
  const auto &matrix_b = GetInput()[1];

  SparseMatrix matrix_b_t = TransposeMatrix(matrix_b);
  int rows = matrix_a.num_rows;

  std::vector<std::vector<double>> temp_values(rows);
  std::vector<std::vector<int>> temp_cols(rows);

#pragma omp parallel for schedule(dynamic, 10) default(none) shared(matrix_a, matrix_b_t, temp_values, temp_cols, rows)
  for (int i = 0; i < rows; ++i) {
    std::vector<double> local_vals;
    std::vector<int> local_cols;

    for (int j = 0; j < matrix_b_t.num_rows; ++j) {
      double sum = DolovVCrsMatMultOmp::DotProduct(matrix_a, i, matrix_b_t, j);
      if (std::fabs(sum) > 1e-15) {
        local_vals.push_back(sum);
        local_cols.push_back(j);
      }
    }
    temp_values[i] = std::move(local_vals);
    temp_cols[i] = std::move(local_cols);
  }

  SparseMatrix res;
  res.num_rows = rows;
  res.num_cols = matrix_b.num_cols;
  res.row_pointers.assign(rows + 1, 0);

  for (int i = 0; i < rows; ++i) {
    res.row_pointers[i + 1] = res.row_pointers[i] + static_cast<int>(temp_values[i].size());
  }

  int total_nz = res.row_pointers[rows];
  res.values.reserve(total_nz);
  res.col_indices.reserve(total_nz);

  for (int i = 0; i < rows; ++i) {
    res.values.insert(res.values.end(), temp_values[i].begin(), temp_values[i].end());
    res.col_indices.insert(res.col_indices.end(), temp_cols[i].begin(), temp_cols[i].end());
  }

  GetOutput() = std::move(res);
  return true;
}

bool DolovVCrsMatMultOmp::PostProcessingImpl() {
  return true;
}

}  // namespace dolov_v_crs_mat_mult
