#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace zyuzin_n_multi_integrals_simpson {

using IntegrandFunc = std::function<double(const std::vector<double> &)>;

struct SimpsonInput {
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  std::vector<int> n_steps;
  IntegrandFunc func;
};

using InType = SimpsonInput;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zyuzin_n_multi_integrals_simpson
