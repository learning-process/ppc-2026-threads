#pragma once

#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<InType, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shkryleva_s_shell_sort_simple_merge
