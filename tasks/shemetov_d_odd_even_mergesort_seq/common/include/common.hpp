#pragma once

#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace shemetov_d_odd_even_mergesort {

using InType = std::tuple<int, std::vector<int>>;
using OutType = std::vector<int>;
using TestType = std::tuple<int, std::vector<int>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shemetov_d_odd_even_mergesort
