#pragma once

#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace zhurin_i_gauss_kernel_seq {

using InType = std::tuple<int, int, int, std::vector<std::vector<int>>>;
// width, height, parts, image

using OutType = std::vector<std::vector<int>>;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

// inline const ppc::task::TaskSettings PPC_SETTINGS_zhurin_i_gauss_kernel_seq{};

}  // namespace zhurin_i_gauss_kernel_seq
