#pragma once

#include <complex>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace barkalova_m_mult_matrix_ccs {

using Complex = std::complex<double>;

struct CCSMatrix {
  std::vector<Complex> values;
  std::vector<int> row_indices;
  std::vector<int> col_ptrs;
  int rows = 0;
  int cols = 0;
  int nnz = 0;
};

using InType = std::pair<CCSMatrix, CCSMatrix>;
using OutType = CCSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;
}  // namespace barkalova_m_mult_matrix_ccs
