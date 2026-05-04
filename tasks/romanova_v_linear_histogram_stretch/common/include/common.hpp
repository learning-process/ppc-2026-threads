#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<bool, size_t, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace romanova_v_linear_histogram_stretch_threads
