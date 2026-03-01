#pragma once

#include <cmath>
#include <vector>

#include "task/include/task.hpp"

namespace orehov_n_jarvis_pass_seq {

struct Point {
  double x;
  double y;
  explicit Point(double x = 0, double y = 0) : x(x), y(y) {}
  bool operator==(const Point &r) const {
    return ((std::abs(x - r.x) < 0.0001) && (std::abs(y - r.y) < 0.0001));
  }
  bool operator<(const Point &other) const {
    if (x != other.x) {
      return x < other.x;
    }
    return y < other.y;
  }
};

using InType = std::vector<Point>;
using OutType = std::vector<Point>;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace orehov_n_jarvis_pass_seq
