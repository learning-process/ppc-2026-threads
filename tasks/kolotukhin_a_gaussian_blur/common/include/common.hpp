#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kolotukhin_a_gaussian_blur {

// pixel data, image wigth, image height
using InType = std::tuple<std::vector<std::uint8_t>, int, int>;
using OutType = std::vector<std::uint8_t>;
using TestType = std::tuple<InType, std::vector<std::uint8_t>, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kolotukhin_a_gaussian_blur
