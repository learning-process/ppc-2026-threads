#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace agafonov_i_matrix_ccs_seq {

struct CCSMatrix {
  size_t rows_num = 0;
  size_t cols_num = 0;
  size_t nnz = 0;
  std::vector<double> vals;
  std::vector<int> row_inds;
  std::vector<int> col_ptrs;
};

using InType = std::pair<CCSMatrix, CCSMatrix>;
using OutType = CCSMatrix;
using BaseTask = ppc::task::Task<InType, OutType>;
using TestType = std::tuple<int, std::string>;

}  // namespace agafonov_i_matrix_ccs_seq
