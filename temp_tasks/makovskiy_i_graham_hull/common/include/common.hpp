#pragma once

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace makovskiy_i_graham_hull {

struct Point {
  double x = 0.0;
  double y = 0.0;

  bool operator==(const Point &other) const {
    return std::abs(x - other.x) < 1e-9 && std::abs(y - other.y) < 1e-9;
  }
};

using InType = std::vector<Point>;
using OutType = std::vector<Point>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace makovskiy_i_graham_hull
