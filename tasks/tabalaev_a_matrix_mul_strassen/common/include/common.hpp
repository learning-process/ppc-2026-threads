#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <cstddef>

#include "task/include/task.hpp"

namespace tabalaev_a_matrix_mul_strassen {

struct MatrixData {
  size_t a_rows{0};
  size_t a_cols_b_rows{0};
  size_t b_cols{0};
  std::vector<double> a;
  std::vector<double> b;
};

using InType = MatrixData;
using OutType = std::vector<double>;
using TestType = std::tuple<size_t, size_t, size_t, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace tabalaev_a_matrix_mul_strassen
