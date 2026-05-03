#pragma once

#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace rozenberg_a_quicksort_simple_merge {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace rozenberg_a_quicksort_simple_merge
