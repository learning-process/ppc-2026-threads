#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace golovanov_d_radix_merge {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::vector<double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace golovanov_d_radix_merge
