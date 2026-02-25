#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
