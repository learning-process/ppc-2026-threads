#pragma once

#include <complex>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

using Complex = std::complex<double>;

struct SparseMatrixCCS {
  int rows{};
  int cols{};

  std::vector<Complex> values{};
  std::vector<int> row_ind{};
  std::vector<int> col_ptr{};
};

using InType = std::tuple<SparseMatrixCCS, SparseMatrixCCS>;
using OutType = SparseMatrixCCS;
using TestType = std::tuple<int, std::string>;

using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
