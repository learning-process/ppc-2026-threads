#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace leonova_a_radix_merge_sort {

using InType = std::vector<int64_t>;
using OutType = std::vector<int64_t>;
using TestType = std::tuple<std::vector<int64_t>, std::vector<int64_t>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace leonova_a_radix_merge_sort
