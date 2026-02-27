#pragma once

#include <complex>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace borunov_v_complex_ccs_seq {

struct SparseMatrix {
  int num_rows = 0;
  int num_cols = 0;
  std::vector<std::complex<double>> values;
  std::vector<int> row_indices;
  std::vector<int> col_ptrs;

  bool operator==(const SparseMatrix &other) const {
    return num_rows == other.num_rows && num_cols == other.num_cols && values == other.values &&
           row_indices == other.row_indices && col_ptrs == other.col_ptrs;
  }
};

using InType = std::vector<SparseMatrix>;
using OutType = std::vector<SparseMatrix>;
using TestType = std::tuple<int, int, int>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace borunov_v_complex_ccs_seq
