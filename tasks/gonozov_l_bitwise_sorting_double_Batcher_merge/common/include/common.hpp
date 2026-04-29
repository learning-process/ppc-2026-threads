#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

using InType = std::vector<double>;   // неотстортированный массив
using OutType = std::vector<double>;  // отстортированный массив
using TestType = std::tuple<std::vector<double>, std::vector<double>, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
