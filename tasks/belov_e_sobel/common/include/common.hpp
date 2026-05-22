#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace belov_e_sobel {

using InType = std::tuple<std::vector<uint8_t>, int, int>;
using OutType = std::tuple<std::vector<uint8_t>, int, int>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace belov_e_sobel
