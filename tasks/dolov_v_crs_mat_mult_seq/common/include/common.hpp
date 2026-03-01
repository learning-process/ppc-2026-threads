#pragma once

#include <string>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace dolov_v_crs_mat_mult_seq {

struct SparseMatrix {
  int num_rows = 0;
  int num_cols = 0;
  std::vector<double> values;
  std::vector<int> col_indices;
  std::vector<int> row_pointers;
};

using InType = std::vector<SparseMatrix>;
using OutType = SparseMatrix;
using TestType = std::pair<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dolov_v_crs_mat_mult_seq
