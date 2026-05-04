#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace tsarkov_k_jarvis_convex_hull {

struct Point {
  int x = 0;
  int y = 0;

  bool operator==(const Point &other_point) const {
    return (x == other_point.x) && (y == other_point.y);
  }
};

using InType = std::vector<Point>;
using OutType = std::vector<Point>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace tsarkov_k_jarvis_convex_hull
