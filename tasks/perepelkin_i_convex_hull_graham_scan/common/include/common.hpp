#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

using InType = std::vector<std::pair<double, double>>;
using OutType = std::vector<std::pair<double, double>>;
using TestType = std::tuple<std::string, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace perepelkin_i_convex_hull_graham_scan
