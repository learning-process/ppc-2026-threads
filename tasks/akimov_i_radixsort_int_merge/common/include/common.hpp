#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace akimov_i_radixsort_int_merge {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace akimov_i_radixsort_int_merge
