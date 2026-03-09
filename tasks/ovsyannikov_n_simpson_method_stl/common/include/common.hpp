#pragma once
#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace ovsyannikov_n_simpson_method_stl {
struct IntegralParams {
  double ax, bx, ay, by;
  int nx, ny;
};
using InType = IntegralParams;
using OutType = double;
using TestType = std::tuple<InType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;
}  // namespace ovsyannikov_n_simpson_method_stl
