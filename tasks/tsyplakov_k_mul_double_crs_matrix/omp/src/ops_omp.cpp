#include "tsyplakov_k_mul_double_crs_matrix/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovKTestTaskOMP::TsyplakovKTestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskOMP::ValidationImpl() {
  const auto &input = GetInput();
  return input.a.cols == input.b.rows;
}

bool TsyplakovKTestTaskOMP::PreProcessingImpl() {
  return true;
}

bool TsyplakovKTestTaskOMP::RunImpl() {
  const auto &input = GetInput();
  const auto &a = input.a;
  const auto &b = input.b;

  const int rows = a.rows;

  std::vector<std::vector<double>> row_values(rows);
  std::vector<std::vector<int>> row_cols(rows);

#pragma omp parallel for schedule(dynamic) default(none) shared(a, b, row_values, row_cols, rows)
  for (int i = 0; i < rows; ++i) {
    std::unordered_map<int, double> acc;

    for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
      int k = a.col_index[j];
      double val_a = a.values[j];

      for (int zz = b.row_ptr[k]; zz < b.row_ptr[k + 1]; ++zz) {
        int j1 = b.col_index[zz];
        acc[j1] += val_a * b.values[zz];
      }
    }

    row_values[i].reserve(acc.size());
    row_cols[i].reserve(acc.size());

    for (const auto &[col, val] : acc) {
      if (std::fabs(val) > 1e-12) {
        row_cols[i].push_back(col);
        row_values[i].push_back(val);
      }
    }
  }

  SparseMatrixCRS c(a.rows, b.cols);

  for (int i = 0; i < c.rows; ++i) {
    c.row_ptr[i + 1] = c.row_ptr[i] + static_cast<int>(row_values[i].size());
  }

  const int nnz = c.row_ptr[c.rows];

  c.values.reserve(nnz);
  c.col_index.reserve(nnz);

  for (int i = 0; i < c.rows; ++i) {
    c.values.insert(c.values.end(), row_values[i].begin(), row_values[i].end());
    c.col_index.insert(c.col_index.end(), row_cols[i].begin(), row_cols[i].end());
  }

  GetOutput() = std::move(c);

  return true;
}

bool TsyplakovKTestTaskOMP::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
