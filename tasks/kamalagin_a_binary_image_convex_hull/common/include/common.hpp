#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kamalagin_a_binary_image_convex_hull {

/// Pixel coordinates for convex hull vertices (x = column, y = row).
struct Point {
  int x = 0;
  int y = 0;

  bool operator==(const Point &other) const {
    return x == other.x && y == other.y;
  }
};

/// Binary image: 0 = background, non-zero = foreground.
struct BinaryImage {
  int rows = 0;
  int cols = 0;
  std::vector<uint8_t> data;

  [[nodiscard]] bool At(int r, int c) const {
    return r >= 0 && r < rows && c >= 0 && c < cols &&
           data[(static_cast<size_t>(r) * static_cast<size_t>(cols)) + static_cast<size_t>(c)] != 0;
  }

  [[nodiscard]] size_t Index(int r, int c) const {
    return (static_cast<size_t>(r) * static_cast<size_t>(cols)) + static_cast<size_t>(c);
  }
};

/// One convex hull = ordered list of vertices (e.g. clockwise).
using Hull = std::vector<Point>;
/// Output: one hull per connected component (same order as component discovery).
using HullList = std::vector<Hull>;

using InType = BinaryImage;
using OutType = HullList;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kamalagin_a_binary_image_convex_hull
