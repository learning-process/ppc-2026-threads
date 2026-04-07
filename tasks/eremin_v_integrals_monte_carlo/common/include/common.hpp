#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace eremin_v_integrals_monte_carlo {

using FunctionType = double (*)(const std::vector<double> &);

struct MonteCarloInput {
  std::vector<std::pair<double, double>> bounds;
  int samples = 0;
  FunctionType func = nullptr;
};

using InType = MonteCarloInput;
using OutType = double;
using TestType =
    std::tuple<int, std::vector<std::pair<double, double>>, int, double (*)(const std::vector<double> &), double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace eremin_v_integrals_monte_carlo
