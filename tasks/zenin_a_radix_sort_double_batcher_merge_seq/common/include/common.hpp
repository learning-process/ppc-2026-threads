#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::tuple<InType, OutType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
