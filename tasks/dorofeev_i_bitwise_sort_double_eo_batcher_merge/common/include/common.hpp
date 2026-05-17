#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

// Изменили int на std::vector<double>
using InType = std::vector<double>;
using OutType = std::vector<double>;

using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
