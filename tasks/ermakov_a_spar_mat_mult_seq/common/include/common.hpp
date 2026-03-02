#pragma once

#include <complex>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace ermakov_a_spar_mat_mult_seq {

struct MatrixCRS {
  int rows = 0;
  int cols = 0;
  std::vector<std::complex<double>> values;
  std::vector<int> col_index;
  std::vector<int> row_ptr;
};

struct InData {
  MatrixCRS A;
  MatrixCRS B;
};

using OutData = MatrixCRS;

using InType = InData;
using OutType = OutData;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace ermakov_a_spar_mat_mult_seq
