#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace sosnina_a_radix_simple_merge {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<size_t, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sosnina_a_radix_simple_merge
