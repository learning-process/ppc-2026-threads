#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace dergynov_s_integrals_multistep_rectangle {

using InType = std::tuple<std::function<double(const std::vector<double>&)>, std::vector<std::pair<double, double>>, int>;
using OutType = double;
using TestType = std::tuple<std::string, InType, double>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dergynov_s_integrals_multistep_rectangle