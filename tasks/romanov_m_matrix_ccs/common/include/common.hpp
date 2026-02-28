#pragma once

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace romanov_m_matrix_ccs {

struct MatrixCCS {
  size_t rows_num = 0;
  size_t cols_num = 0;
  size_t nnz = 0;
  std::vector<double> vals;
  std::vector<size_t> row_inds;
  std::vector<size_t> col_ptrs;
};

using InType = std::pair<MatrixCCS, MatrixCCS>;
using OutType = MatrixCCS;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace romanov_m_matrix_ccs
