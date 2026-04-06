#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs_seq {

struct CCSMatrix {
  std::vector<double> values;
  std::vector<size_t> row_indices;
  std::vector<size_t> col_ptrs;
  size_t rows = 0;
  size_t cols = 0;
  size_t nnz = 0;
};

using InType = std::pair<CCSMatrix, CCSMatrix>;
using OutType = CCSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kapanova_s_sparse_matrix_mult_ccs_seq
