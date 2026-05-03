#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace egorova_l_binary_convex_hull {

struct Point {
  int x, y;
  bool operator<(const Point &other) const {
    return std::tie(x, y) < std::tie(other.x, other.y);
  }
};

struct ImageData {
  std::vector<uint8_t> data;
  int width = 0;
  int height = 0;
};

using InType = ImageData;
using OutType = std::vector<std::vector<Point>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace egorova_l_binary_convex_hull
