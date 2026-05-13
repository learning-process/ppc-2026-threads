#pragma once

#include <complex>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_all {

/// Compressed columns (CCS): for each j in [0, col_count), nonzeros in column j are
/// (row_index[t], nz[t]) with t in [col_start[j], col_start[j+1]) and row_index strictly increasing.
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

}  // namespace sabutay_sparse_complex_ccs_mult_all
