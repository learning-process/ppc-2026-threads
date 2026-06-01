#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace gusev_d_double_sort_even_odd_batcher_tbb_task_threads {

using ValueType = double;
using InType = std::vector<ValueType>;
using OutType = std::vector<ValueType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace gusev_d_double_sort_even_odd_batcher_tbb_task_threads
