#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace muhammadkhon_i_stressen_alg {

struct MatrixInput {
  std::vector<double> a;
  std::vector<double> b;
  int n{};
};

using InType = MatrixInput;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace muhammadkhon_i_stressen_alg
