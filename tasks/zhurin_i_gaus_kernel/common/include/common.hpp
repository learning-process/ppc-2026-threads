#pragma once

#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace zhurin_i_gaus_kernel {

using InType = std::tuple<int, int, int, std::vector<std::vector<int>>>;
using OutType = std::vector<std::vector<int>>;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

// inline const ppc::task::TaskSettings PPC_SETTINGS_zhurin_i_gaus_kernel{};

}  // namespace zhurin_i_gaus_kernel
