#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace badanov_a_select_edge_sobel_seq {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace badanov_a_select_edge_sobel_seq
