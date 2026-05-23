#pragma once
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace solonin_v_radix_sort_batcher {
using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;
}  // namespace solonin_v_radix_sort_batcher
