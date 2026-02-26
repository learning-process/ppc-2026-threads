#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace golovanov_d_radix_merge {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace golovanov_d_radix_merge
