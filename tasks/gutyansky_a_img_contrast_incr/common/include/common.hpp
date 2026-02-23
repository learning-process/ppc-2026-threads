#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace gutyansky_a_img_contrast_incr {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace gutyansky_a_img_contrast_incr
