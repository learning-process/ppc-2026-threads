#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

struct SparseMatrixCRS {
  int rows;
  int cols;
  std::vector<double> values;
  std::vector<int> col_index;
  std::vector<int> row_ptr;

  SparseMatrixCRS() : rows(0), cols(0) {}

  SparseMatrixCRS(int n_rows, int n_cols) : rows(n_rows), cols(n_cols), row_ptr(n_rows + 1, 0) {}
};

struct MatrixMulInput {
  SparseMatrixCRS A;
  SparseMatrixCRS B;
};

using InType = MatrixMulInput;
using OutType = SparseMatrixCRS;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace tsyplakov_k_mul_double_crs_matrix
