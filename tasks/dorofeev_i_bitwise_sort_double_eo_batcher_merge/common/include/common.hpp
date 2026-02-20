#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
