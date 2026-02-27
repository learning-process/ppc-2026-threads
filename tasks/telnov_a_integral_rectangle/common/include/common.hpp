#pragma once

#include <string>
#include <tuple>
#include <utility>

namespace telnov_a_integral_rectangle {

using InType = std::pair<int, int>;
using OutType = double;
using TestType = std::tuple<InType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace telnov_a_integral_rectangle
