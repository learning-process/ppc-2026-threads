#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sokolov_k_matrix_double_fox_seq {

using InType = std::tuple<int, int, std::vector<double>, std::vector<double>>;
using OutType = std::vector<double>;
using TestType = std::tuple<std::string, int, int, std::vector<double>, std::vector<double>, std::vector<double>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sokolov_k_matrix_double_fox_seq
