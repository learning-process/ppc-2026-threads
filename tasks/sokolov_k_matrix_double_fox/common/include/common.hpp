#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace sokolov_k_matrix_double_fox {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sokolov_k_matrix_double_fox
