#pragma once

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {

struct SparseMatrixCCS {
  std::vector<double> values;
  std::vector<int> row_indices;
  std::vector<int> col_ptrs;
  int rows;
  int cols;

  SparseMatrixCCS() : rows(0), cols(0) {}

  SparseMatrixCCS(int n_rows, int n_cols) : col_ptrs(n_cols + 1, 0), rows(n_rows), cols(n_cols) {}
};

using InType = std::pair<SparseMatrixCCS, SparseMatrixCCS>;
using OutType = SparseMatrixCCS;
using TestType = std::tuple<int, int, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
