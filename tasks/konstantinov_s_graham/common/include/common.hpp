#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace konstantinov_a_graham {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace konstantinov_a_graham
