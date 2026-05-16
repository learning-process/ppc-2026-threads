#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sinev_a_mult_matrix_fox_algorithm_seq {

using InType = std::tuple<size_t, std::vector<double>, std::vector<double>>;
using OutType = std::vector<double>;
using TestType = std::tuple<size_t, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sinev_a_mult_matrix_fox_algorithm_seq
