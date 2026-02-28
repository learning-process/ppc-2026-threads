#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace litvyakov_d_shell_sort {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace litvyakov_d_shell_sort
