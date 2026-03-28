#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge_seq {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<InType, OutType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace ovchinnikov_m_shell_sort_batcher_merge_seq
