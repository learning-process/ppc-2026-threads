#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace korolev_k_matrix_mult {

struct MatrixInput {
  size_t n;
  std::vector<double> A;  // row-major, n*n
  std::vector<double> B;
};

using InType = MatrixInput;
using OutType = std::vector<double>;  // C = A*B, row-major, n*n
using TestType = std::tuple<size_t, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace korolev_k_matrix_mult
