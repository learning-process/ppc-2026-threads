#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace papulina_y_radix_sort {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::tuple<std::vector<double>, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace papulina_y_radix_sort
