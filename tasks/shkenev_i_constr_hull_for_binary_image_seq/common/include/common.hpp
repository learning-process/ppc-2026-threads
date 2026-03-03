#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace shkenev_i_constr_hull_for_binary_image_seq {

struct Point {
  int x{0};
  int y{0};

  Point() = default;
  Point(int px, int py) : x(px), y(py) {}

  bool operator==(const Point &other) const {
    return x == other.x && y == other.y;
  }

  bool operator!=(const Point &other) const {
    return !(*this == other);
  }

  bool operator<(const Point &other) const {
    return (x < other.x) || (x == other.x && y < other.y);
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

}  // namespace shkenev_i_constr_hull_for_binary_image_seq
