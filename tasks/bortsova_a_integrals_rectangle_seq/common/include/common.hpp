#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace bortsova_a_integrals_rectangle_seq {

using Function = std::function<double(const std::vector<double> &)>;

struct IntegralInput {
  Function func;
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  int num_steps = 0;
};

using InType = IntegralInput;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace bortsova_a_integrals_rectangle_seq
