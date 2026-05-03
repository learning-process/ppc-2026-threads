#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

struct SparseMatCCS {
  int rows_count = 0;
  int cols_count = 0;
  int non_zeros = 0;
  std::vector<int> col_ptrs;
  std::vector<int> row_indices;
  std::vector<double> values;

  SparseMatCCS() = default;
  SparseMatCCS(const SparseMatCCS &) = default;
  SparseMatCCS &operator=(const SparseMatCCS &) = default;
};

using InType = std::tuple<SparseMatCCS, SparseMatCCS>;
using OutType = SparseMatCCS;
using TestType = std::tuple<std::string, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace volkov_a_sparse_mat_mul_ccs
