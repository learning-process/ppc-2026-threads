#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

using InType = std::pair<std::vector<double>, std::vector<double>>;
using OutType = std::vector<double>;
using TestType = std::tuple<InType, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
