#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace samoylenko_i_integral_trapezoid {

struct Input {
  std::vector<double> a;
  std::vector<double> b;
  std::vector<int> n;
  int function_choice{};
};

using InType = Input;
using OutType = double;
using TestType = std::pair<InType, double>;

using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace samoylenko_i_integral_trapezoid
