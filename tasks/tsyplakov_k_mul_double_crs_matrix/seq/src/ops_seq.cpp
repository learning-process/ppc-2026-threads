#include "tsyplakov_k_mul_double_crs_matrix/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <vector>

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovKTestTaskSEQ::TsyplakovKTestTaskSEQ(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskSEQ::ValidationImpl() {
  return (GetInput().a.cols == GetInput().b.rows) && (GetInput().a.rows > 0) && (GetInput().a.cols > 0) &&
         (GetInput().b.rows > 0) && (GetInput().b.cols > 0);
}

bool TsyplakovKTestTaskSEQ::PreProcessingImpl() {
  GetOutput() = SparseMatrixCRS(GetInput().a.rows, GetInput().b.cols);
  return true;
}

std::vector<double> TsyplakovKTestTaskSEQ::MultiplyRowByMatrix(const std::vector<double>& row_values,
                                                               const std::vector<int>& row_cols,
                                                               const SparseMatrixCRS& b, int& result_nnz) {
  std::map<int, double> temp_result;
  for (size_t i = 0; i < row_cols.size(); ++i) {
    int col_a = row_cols[i];
    double val_a = row_values[i];

    int start = b.row_ptr[col_a];
    int end = b.row_ptr[col_a + 1];

    for (int j = start; j < end; ++j) {
      int col_b = b.col_index[j];
      double val_b = b.values[j];

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
  const SparseMatrixCRS& a = GetInput().a;
  const SparseMatrixCRS& b = GetInput().b;
  SparseMatrixCRS& c = GetOutput();

  c.values.clear();
  c.col_index.clear();
  c.row_ptr.assign(a.rows + 1, 0);

  for (int i = 0; i < a.rows; ++i) {
    int start_a = a.row_ptr[i];
    int end_a = a.row_ptr[i + 1];

    std::vector<double> row_values(a.values.begin() + start_a, a.values.begin() + end_a);
    std::vector<int> row_cols(a.col_index.begin() + start_a, a.col_index.begin() + end_a);

    int row_nnz = 0;
    std::vector<double> row_result = MultiplyRowByMatrix(row_values, row_cols, b, row_nnz);

    c.row_ptr[i + 1] = c.row_ptr[i] + row_nnz;
  }

  c.values.clear();
  c.col_index.clear();
  c.row_ptr.assign(a.rows + 1, 0);

  for (int i = 0; i < a.rows; ++i) {
    int start_a = a.row_ptr[i];
    int end_a = a.row_ptr[i + 1];

    std::map<int, double> temp_result;

    for (int k_idx = start_a; k_idx < end_a; ++k_idx) {
      int k = a.col_index[k_idx];
      double a_ik = a.values[k_idx];

      int start_b = b.row_ptr[k];
      int end_b = b.row_ptr[k + 1];

      for (int j_idx = start_b; j_idx < end_b; ++j_idx) {
        int j = b.col_index[j_idx];
        double b_kj = b.values[j_idx];
        temp_result[j] += a_ik * b_kj;
      }
    }

    const double eps = std::numeric_limits<double>::epsilon() * 100;
    c.row_ptr[i] = static_cast<int>(c.values.size());

    for (const auto& [col, val] : temp_result) {
      if (std::abs(val) > eps) {
        c.col_index.push_back(col);
        c.values.push_back(val);
      }
    }
  }
  c.row_ptr[a.rows] = static_cast<int>(c.values.size());

  return true;
}

bool TsyplakovKTestTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
