#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace yurkin_g_graham_scan {

struct Point {
  double x;
  double y;
  bool operator==(const Point &other) const {
    return x == other.x && y == other.y;
  }
};

using InType = std::vector<Point>;
using OutType = std::vector<Point>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace yurkin_g_graham_scan
