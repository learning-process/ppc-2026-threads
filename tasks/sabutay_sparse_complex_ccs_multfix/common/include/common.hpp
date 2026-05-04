#pragma once

#include <complex>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_multfix {
struct CCS {
  int row_count{0};
  int col_count{0};
  std::vector<int> col_start;
  std::vector<int> row_index;
  std::vector<std::complex<double>> nz;

  CCS() = default;
  CCS(const CCS &) = default;
  CCS &operator=(const CCS &) = default;
};

using InType = std::tuple<CCS, CCS>;
using OutType = CCS;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sabutay_sparse_complex_ccs_multfix
