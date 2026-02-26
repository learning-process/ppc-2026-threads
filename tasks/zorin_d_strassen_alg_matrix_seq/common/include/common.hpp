#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

struct MatMulInput {
  std::size_t n{};                 
  std::vector<double> A;           
  std::vector<double> B;           
};

using InType = MatMulInput;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zorin_d_strassen_alg_matrix_seq
