#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace nikolaev_d_block_linear_image_filtering {

using InType = std::tuple<int, int, std::vector<uint8_t>>;  // img width, img height, pxl data
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<InType, OutType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nikolaev_d_block_linear_image_filtering
