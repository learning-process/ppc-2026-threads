#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sizov_d_sparse_crs_mult {

struct CRSMatrix {
  std::vector<double> values;
  std::vector<std::size_t> col_indices;
  std::vector<std::size_t> row_ptr;
  std::size_t rows = 0;
  std::size_t cols = 0;
};

using InType = std::tuple<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using TestType = std::tuple<std::string, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
                            std::vector<std::vector<double>>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sizov_d_sparse_crs_mult
