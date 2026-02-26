#include "tsyplakov_k_mul_double_crs_matrix/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovKTestTaskSEQ::TsyplakovKTestTaskSEQ(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskSEQ::ValidationImpl() {
  return (GetInput().A.cols == GetInput().B.rows) && (GetInput().A.rows > 0) && (GetInput().A.cols > 0) &&
         (GetInput().B.rows > 0) && (GetInput().B.cols > 0);
}

bool TsyplakovKTestTaskSEQ::PreProcessingImpl() {
  GetOutput() = SparseMatrixCRS(GetInput().A.rows, GetInput().B.cols);
  return true;
}

std::vector<double> TsyplakovKTestTaskSEQ::MultiplyRowByMatrix(const std::vector<double>& row_values,
                                                               const std::vector<int>& row_cols,
                                                               const SparseMatrixCRS& B, int& result_nnz) {
  std::map<int, double> temp_result;
  for (size_t i = 0; i < row_cols.size(); ++i) {
    int col_a = row_cols[i];
    double val_a = row_values[i];

    int start = B.row_ptr[col_a];
    int end = B.row_ptr[col_a + 1];

    for (int j = start; j < end; ++j) {
      int col_b = B.col_index[j];
      double val_b = B.values[j];

      temp_result[col_b] += val_a * val_b;
    }
  }

  const double eps = std::numeric_limits<double>::epsilon() * 100;
  std::vector<double> result_values;
  result_nnz = 0;

  for (auto it = temp_result.begin(); it != temp_result.end();) {
    if (std::abs(it->second) < eps) {
      it = temp_result.erase(it);
    } else {
      result_values.push_back(it->second);
      ++it;
      ++result_nnz;
    }
  }

  return result_values;
}

bool TsyplakovKTestTaskSEQ::RunImpl() {
  const SparseMatrixCRS& A = GetInput().A;
  const SparseMatrixCRS& B = GetInput().B;
  SparseMatrixCRS& C = GetOutput();

  C.values.clear();
  C.col_index.clear();
  C.row_ptr.assign(A.rows + 1, 0);

  for (int i = 0; i < A.rows; ++i) {
    int start_a = A.row_ptr[i];
    int end_a = A.row_ptr[i + 1];

    std::vector<double> row_values(A.values.begin() + start_a, A.values.begin() + end_a);
    std::vector<int> row_cols(A.col_index.begin() + start_a, A.col_index.begin() + end_a);

    int row_nnz = 0;
    std::vector<double> row_result = MultiplyRowByMatrix(row_values, row_cols, B, row_nnz);

    C.row_ptr[i + 1] = C.row_ptr[i] + row_nnz;
  }

  C.values.clear();
  C.col_index.clear();
  C.row_ptr.assign(A.rows + 1, 0);

  for (int i = 0; i < A.rows; ++i) {
    int start_a = A.row_ptr[i];
    int end_a = A.row_ptr[i + 1];

    std::map<int, double> temp_result;

    for (int k_idx = start_a; k_idx < end_a; ++k_idx) {
      int k = A.col_index[k_idx];
      double a_ik = A.values[k_idx];

      int start_b = B.row_ptr[k];
      int end_b = B.row_ptr[k + 1];

      for (int j_idx = start_b; j_idx < end_b; ++j_idx) {
        int j = B.col_index[j_idx];
        double b_kj = B.values[j_idx];
        temp_result[j] += a_ik * b_kj;
      }
    }

    const double eps = std::numeric_limits<double>::epsilon() * 100;
    C.row_ptr[i] = static_cast<int>(C.values.size());

    for (const auto& [col, val] : temp_result) {
      if (std::abs(val) > eps) {
        C.col_index.push_back(col);
        C.values.push_back(val);
      }
    }
  }
  C.row_ptr[A.rows] = static_cast<int>(C.values.size());

  return true;
}

bool TsyplakovKTestTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
