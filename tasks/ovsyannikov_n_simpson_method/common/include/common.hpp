#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace ovsyannikov_n_simpson_method {

struct IntegralParams {
  double ax, bx;  // Грaницы по X
  double ay, by;  // Грaницы по Y
  int nx, ny;
};

using InType = IntegralParams;
using OutType = double;
using TestType = std::tuple<InType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace ovsyannikov_n_simpson_method
