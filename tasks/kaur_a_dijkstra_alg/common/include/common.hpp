#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace kaur_a_dijkstra_alg {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kaur_a_dijkstra_alg
