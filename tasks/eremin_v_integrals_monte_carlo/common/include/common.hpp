#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace eremin_v_integrals_monte_carlo {

using FunctionType = std::function<double(const std::vector<double> &)>;

struct MonteCarloInput {
  std::vector<std::pair<double, double>> bounds;
  int samples;
  FunctionType func;
};

using InType = MonteCarloInput;
using OutType = double;
using TestType = std::tuple<MonteCarloInput, double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace eremin_v_integrals_monte_carlo
