#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace dilshodov_a_spmm_double_css {

using DenseMatrix = std::vector<std::vector<double>>;

struct SparseMatrixCCS {
  int rows_count = 0;
  int cols_count = 0;
  int non_zeros = 0;
  std::vector<int> col_ptrs;
  std::vector<int> row_indices;
  std::vector<double> values;
};

using InType = std::tuple<SparseMatrixCCS, SparseMatrixCCS>;
using OutType = SparseMatrixCCS;
using TestType = std::tuple<std::string, DenseMatrix, DenseMatrix>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dilshodov_a_spmm_double_css
