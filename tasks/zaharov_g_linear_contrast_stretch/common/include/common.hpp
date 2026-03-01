#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace zaharov_g_linear_contrast_stretch {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zaharov_g_linear_contrast_stretch
