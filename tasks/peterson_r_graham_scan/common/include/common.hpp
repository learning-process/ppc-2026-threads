#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace peterson_r_graham_scan {

using InputValue = int;
using OutputValue = int;
using TestParameters = std::tuple<int, std::string>;
using TaskBase = ppc::task::Task<InputValue, OutputValue>;

struct Point2D {
  double coord_x;
  double coord_y;

  Point2D() = default;
  Point2D(double x_val, double y_val) : coord_x(x_val), coord_y(y_val) {}
};

using PointSet = std::vector<Point2D>;

}  // namespace peterson_r_graham_scan
