#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace otcheskov_s_contrast_lin_stretch {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<std::string, size_t>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace otcheskov_s_contrast_lin_stretch
