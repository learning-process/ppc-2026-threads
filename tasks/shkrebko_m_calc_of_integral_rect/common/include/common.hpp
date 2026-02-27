#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace shkrebko_m_calc_of_integral_rect {

struct IntegralInput {
  std::vector<std::pair<double, double>> limits;
  int steps = 0;
  std::function<double(const std::vector<double> &)> func;
};

using InType = IntegralInput;
using OutType = double;
using TestType = std::tuple<std::string, InType, double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shkrebko_m_calc_of_integral_rect
