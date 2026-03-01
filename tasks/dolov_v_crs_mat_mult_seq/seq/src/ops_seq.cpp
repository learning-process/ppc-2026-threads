#include "dolov_v_crs_mat_mult_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "dolov_v_crs_mat_mult_seq/common/include/common.hpp"

namespace dolov_v_crs_mat_mult_seq {

DolovVCrsMatMultSeq::DolovVCrsMatMultSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DolovVCrsMatMultSeq::ValidationImpl() {
  const auto &input_data = GetInput();
  if (input_data.size() != 2) {
    return false;
  }
  const auto &matrix_a = input_data[0];
  const auto &matrix_b = input_data[1];
  return matrix_a.num_cols == matrix_b.num_rows && matrix_a.num_rows > 0 && matrix_b.num_cols > 0;
}

bool DolovVCrsMatMultSeq::PreProcessingImpl() {
  const auto &input_data = GetInput();
  auto &result = GetOutput();
  result.num_rows = input_data[0].num_rows;
  result.num_cols = input_data[1].num_cols;
  result.row_pointers.assign(result.num_rows + 1, 0);
  result.values.clear();
  result.col_indices.clear();
  return true;
}

SparseMatrix DolovVCrsMatMultSeq::TransposeMatrix(const SparseMatrix &matrix) {
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

double DolovVCrsMatMultSeq::DotProduct(const SparseMatrix &matrix_a, int row_a, const SparseMatrix &matrix_b_t,
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

bool DolovVCrsMatMultSeq::RunImpl() {
  const auto &input_data = GetInput();
  const auto &matrix_a = input_data[0];
  const auto &matrix_b = input_data[1];
  auto &result = GetOutput();

  SparseMatrix matrix_b_t = TransposeMatrix(matrix_b);

  for (int i = 0; i < matrix_a.num_rows; ++i) {
    for (int j = 0; j < matrix_b_t.num_rows; ++j) {
      double sum = DotProduct(matrix_a, i, matrix_b_t, j);

      if (std::abs(sum) > 1e-15) {
        result.values.push_back(sum);
        result.col_indices.push_back(j);
      }
    }
    result.row_pointers[i + 1] = static_cast<int>(result.values.size());
  }
  return true;
}

bool DolovVCrsMatMultSeq::PostProcessingImpl() {
  return true;
}

}  // namespace dolov_v_crs_mat_mult_seq
