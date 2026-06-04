#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace nazarova_k_calc_integ_rectangles {

struct IntegrationInput {
  std::function<double(const std::vector<double> &)> function;
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  std::vector<std::size_t> steps;
};

using InType = IntegrationInput;
using OutType = double;
using TestType = std::tuple<InType, OutType, double, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nazarova_k_calc_integ_rectangles
