#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kapanova_s_increase_contrast_seq {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kapanova_s_increase_contrast_seq
