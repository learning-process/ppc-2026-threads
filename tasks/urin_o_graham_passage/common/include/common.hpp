#pragma once

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace urin_o_graham_passage {

struct Point {
  double x;
  double y;

  Point() : x(0.0), y(0.0) {}
  Point(double x_val, double y_val) : x(x_val), y(y_val) {}  // Исправлено: x_, y_ -> x_val, y_val

  bool operator==(const Point &other) const {
    const double eps = 1e-10;
    return std::abs(x - other.x) < eps && std::abs(y - other.y) < eps;
  }

  bool operator!=(const Point &other) const {
    return !(*this == other);
  }
};

using InType = std::vector<Point>;
using OutType = std::vector<Point>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace urin_o_graham_passage
