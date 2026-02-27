#pragma once
#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace buzulukski_d_gaus_gorizontal {
using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;
}  // namespace buzulukski_d_gaus_gorizontal
