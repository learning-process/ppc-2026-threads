#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

using InType = std::tuple<size_t, std::vector<double>, std::vector<double>>;
using OutType = std::vector<double>;
using TestType = std::tuple<size_t, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zyazeva_s_matrix_mult_cannon_alg
