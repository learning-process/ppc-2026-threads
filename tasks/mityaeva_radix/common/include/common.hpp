#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"

namespace mityaeva_radix {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::size_t;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace mityaeva_radix
