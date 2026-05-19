#pragma once

#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace nazyrov_a_multidim_integral_rectangle {

using Func = std::function<double(const std::vector<double> &)>;
using Bounds = std::vector<std::pair<double, double>>;
using InType = std::tuple<Func, Bounds, int>;
using OutType = double;
using TestType = std::tuple<int, int, int>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nazyrov_a_multidim_integral_rectangle
