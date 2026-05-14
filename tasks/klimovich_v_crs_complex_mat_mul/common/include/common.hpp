#pragma once

#include <complex>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace klimovich_v_crs_complex_mat_mul {

using Cplx = std::complex<double>;

constexpr double kZeroDropTol = 1e-12;

struct CrsMatrix {
  int n_rows = 0;
  int n_cols = 0;
  std::vector<int> row_offsets;
  std::vector<int> col_indices;
  std::vector<Cplx> data;

  CrsMatrix() = default;

  CrsMatrix(int rows, int cols) : n_rows(rows), n_cols(cols), row_offsets(rows + 1, 0) {}
};

using InType = std::tuple<CrsMatrix, CrsMatrix>;
using OutType = CrsMatrix;
using TestType = std::tuple<int, int, int>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace klimovich_v_crs_complex_mat_mul
