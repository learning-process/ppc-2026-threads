#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads
