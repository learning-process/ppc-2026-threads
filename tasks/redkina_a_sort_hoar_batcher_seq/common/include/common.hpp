// common.hpp
#pragma once

#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace redkina_a_sort_hoar_batcher_seq {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<int, std::vector<int>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace redkina_a_sort_hoar_batcher_seq
