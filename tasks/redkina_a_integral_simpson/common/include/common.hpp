#pragma once

#include <functional>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace redkina_a_integral_simpson {

struct InputData {
  std::function<double(const std::vector<double> &)> func;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<int> n;
};

using InType = InputData;
using OutType = double;
using TestType = std::tuple<int, double (*)(const std::vector<double> &), std::vector<double>, std::vector<double>,
                            std::vector<int>, double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace redkina_a_integral_simpson
