#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace romanov_a_gauss_block {

using InType = std::tuple<int, int, std::vector<uint8_t>>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<int, int, std::vector<uint8_t>, std::vector<uint8_t>, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace romanov_a_gauss_block
