#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace boltenkov_s_gaussian_kernel {

using InType = std::tuple<std::size_t, std::size_t, std::vector<std::vector<int>>>;
using OutType = std::vector<std::vector<int>>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace boltenkov_s_gaussian_kernel
