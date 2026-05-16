#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

struct InputIntegralData {
  std::vector<double> left_bounds;
  std::vector<double> right_bounds;
  std::vector<int> step_n_size;
  int type_function = 0;
  double epsilon = 1e-2;
};

using InType = InputIntegralData;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
