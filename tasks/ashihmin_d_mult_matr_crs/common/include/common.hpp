#pragma once

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace ashihmin_d_mult_matr_crs {

struct CRSMatrix {
  int rows{};
  int cols{};
  std::vector<double> values;
  std::vector<int> col_index;
  std::vector<int> row_ptr;
};

using InType = std::pair<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;

using DenseMatrix = std::vector<std::vector<double>>;

using TestType = std::tuple<std::string, DenseMatrix, DenseMatrix, DenseMatrix>;

using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace ashihmin_d_mult_matr_crs
