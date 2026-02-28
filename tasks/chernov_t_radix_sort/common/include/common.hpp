#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<std::string, std::vector<int>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace chernov_t_radix_sort
