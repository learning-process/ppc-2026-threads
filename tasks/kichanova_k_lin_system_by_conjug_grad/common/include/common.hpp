#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

struct LinSystemData {
  std::vector<double> A;
  std::vector<double> b;
  std::vector<double> x;
  int n = 0;
  double epsilon = 1e-10;
};

using InType = LinSystemData;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kichanova_k_lin_system_by_conjug_grad
