#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

using InType = int;     // update
using OutType = int;    // update
using TestType = std::tuple<int, std::string>;      // update
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace perepelkin_i_convex_hull_graham_scan
