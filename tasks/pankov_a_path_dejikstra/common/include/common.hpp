#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace pankov_a_path_dejikstra {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace pankov_a_path_dejikstra
