#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // zenin_a_radix_sort_double_batcher_merge_seq
