#pragma once

#include <vector>
#include <string>
#include <tuple>
#include "task/include/task.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging