#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

using ValType = int;
using InType = std::vector<ValType>;
using OutType = std::vector<ValType>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace vasiliev_m_shell_sort_batcher_merge
