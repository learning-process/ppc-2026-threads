#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace rysev_m_linear_filter_gauss_kernel
