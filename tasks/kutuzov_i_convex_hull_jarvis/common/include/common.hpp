#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kutuzov_i_convex_hull_jarvis {

// List of all points
using InType = std::vector<std::tuple<double, double>>;

// Sequence of points that make up the convex hull
using OutType = std::vector<std::tuple<double, double>>;

// List of all points and the expected convex hull
using TestType = std::tuple<InType, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kutuzov_i_convex_hull_jarvis
