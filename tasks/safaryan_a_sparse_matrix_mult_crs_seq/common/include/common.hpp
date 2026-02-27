#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {

struct CRSMatrix {
  std::vector<double> values;
  std::vector<size_t> col_indices;
  std::vector<size_t> row_ptr;
  size_t rows = 0;
  size_t cols = 0;
  size_t nnz = 0;
};

using InType = std::pair<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
