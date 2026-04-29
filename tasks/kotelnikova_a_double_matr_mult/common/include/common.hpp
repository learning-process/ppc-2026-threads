#pragma once

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace kotelnikova_a_double_matr_mult {

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

}  // namespace kotelnikova_a_double_matr_mult
