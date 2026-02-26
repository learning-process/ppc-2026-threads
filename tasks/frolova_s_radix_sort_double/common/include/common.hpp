#pragma once

#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace frolova_s_radix_sort_double