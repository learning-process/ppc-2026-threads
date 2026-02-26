#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace lazareva_a_matrix_mult_strassen {

struct MatrixInput {
  std::vector<double> a;
  std::vector<double> b;
  int n{};
};

using InType = MatrixInput;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace lazareva_a_matrix_mult_strassen