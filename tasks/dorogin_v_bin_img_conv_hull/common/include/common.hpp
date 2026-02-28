#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull {

struct Point {
  int x{0};
  int y{0};

  Point() = default;
  Point(int px, int py) : x(px), y(py) {}

  bool operator==(const Point &oth) const {
    return x == oth.x && y == oth.y;
  }

  bool operator!=(const Point &oth) const {
    return !(*this == oth);
  }

  bool operator<(const Point &oth) const {
    if (x == oth.x) {
      return y < oth.y;
    }
    return x < oth.x;
  }
};

struct BinaryImage {
  int width{0};
  int height{0};
  std::vector<uint8_t> pixels;
  std::vector<std::vector<Point>> components;
  std::vector<std::vector<Point>> convex_hulls;
};

using InType = BinaryImage;
using OutType = BinaryImage;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dorogin_v_bin_img_conv_hull
