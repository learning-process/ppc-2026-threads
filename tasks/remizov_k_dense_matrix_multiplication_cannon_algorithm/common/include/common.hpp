#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

using InType = std::tuple<int, std::vector<std::vector<double>>, std::vector<std::vector<double>>>;
using OutType = std::vector<std::vector<double>>;
using TestType = std::tuple<std::string, int, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
                            std::vector<std::vector<double>>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
