#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace afanasyev_a_integ_rect_method {

using InType = int;
using OutType = double;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace afanasyev_a_integ_rect_method
