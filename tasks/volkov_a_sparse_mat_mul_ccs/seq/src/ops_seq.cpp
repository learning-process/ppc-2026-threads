#include "volkov_a_sparse_mat_mul_ccs/seq/include/ops_seq.hpp"

#include <cmath>
#include <tuple>
#include <vector>

#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

VolkovASparseMatMulCcsSeq::VolkovASparseMatMulCcsSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VolkovASparseMatMulCcsSeq::ValidationImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  return (matrix_a.cols_count == matrix_b.rows_count);
}

bool VolkovASparseMatMulCcsSeq::PreProcessingImpl() {
  return true;
}

bool VolkovASparseMatMulCcsSeq::RunImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  matrix_c.rows_count = matrix_a.rows_count;
  matrix_c.cols_count = matrix_b.cols_count;
  matrix_c.col_ptrs.assign(matrix_c.cols_count + 1, 0);
  matrix_c.row_indices.clear();
  matrix_c.values.clear();

  std::vector<double> col_accumulator(matrix_a.rows_count, 0.0);

  for (int j = 0; j < matrix_b.cols_count; ++j) {
    matrix_c.col_ptrs[j] = static_cast<int>(matrix_c.values.size());

    int b_start = matrix_b.col_ptrs[j];
    int b_end = matrix_b.col_ptrs[j + 1];

    for (int k = b_start; k < b_end; ++k) {
      int b_row = matrix_b.row_indices[k];
      double b_val = matrix_b.values[k];

      int a_start = matrix_a.col_ptrs[b_row];
      int a_end = matrix_a.col_ptrs[b_row + 1];

      for (int idx = a_start; idx < a_end; ++idx) {
        int a_row = matrix_a.row_indices[idx];
        double a_val = matrix_a.values[idx];
        
        col_accumulator[a_row] += a_val * b_val;
      }
    }

    for (int i = 0; i < matrix_a.rows_count; ++i) {
      if (std::abs(col_accumulator[i]) > 1e-10) {
        matrix_c.row_indices.push_back(i);
        matrix_c.values.push_back(col_accumulator[i]);
      }
      col_accumulator[i] = 0.0;
    }
  }

  matrix_c.non_zeros = static_cast<int>(matrix_c.values.size());
  matrix_c.col_ptrs[matrix_b.cols_count] = matrix_c.non_zeros;

  return true;
}

bool VolkovASparseMatMulCcsSeq::PostProcessingImpl() {
  return true;
}

}  // namespace volkov_a_sparse_mat_mul_ccs