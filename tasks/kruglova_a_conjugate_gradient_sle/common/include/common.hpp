#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kruglova_a_conjugate_gradient_sle {

struct SLEInput {
  int size;
  std::vector<double> A;
  std::vector<double> b;
};

using InType = SLEInput;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kruglova_a_conjugate_gradient_sle
