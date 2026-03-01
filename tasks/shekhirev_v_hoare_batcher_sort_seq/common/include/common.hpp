#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shekhirev_v_hoare_batcher_sort_seq
