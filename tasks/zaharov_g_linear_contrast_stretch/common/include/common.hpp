#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace zaharov_g_linear_contrast_stretch {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zaharov_g_linear_contrast_stretch
