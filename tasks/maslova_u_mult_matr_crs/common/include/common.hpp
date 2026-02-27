#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace maslova_u_mult_matr_crs {

struct CRSMatrix {
  int rows = 0;
  int cols = 0;
  std::vector<double> values;
  std::vector<int> col_ind;
  std::vector<int> row_ptr;
};

using InType = std::tuple<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace maslova_u_mult_matr_crs
