#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace chernov_t_radix_sort
