#pragma once

#include <complex>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult {

struct CCS {
  int m = 0;
  int n = 0;
  std::vector<int> col_ptr;
  std::vector<int> row_ind;
  std::vector<std::complex<double>> values;

  CCS() = default;
  CCS(const CCS &) = default;
  CCS &operator=(const CCS &) = default;
};

using InType = std::tuple<CCS, CCS>;
using OutType = CCS;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sabutay_sparse_complex_ccs_mult
