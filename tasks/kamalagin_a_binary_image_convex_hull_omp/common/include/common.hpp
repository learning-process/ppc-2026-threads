#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace kamalagin_a_binary_image_convex_hull_omp {

struct Point {
  int x{};
  int y{};
};

inline bool operator==(const Point &a, const Point &b) {
  return a.x == b.x && a.y == b.y;
}

using Hull = std::vector<Point>;
using HullList = std::vector<Hull>;

struct BinaryImage {
  int rows{};
  int cols{};
  std::vector<uint8_t> data;

  [[nodiscard]] size_t Index(int r, int c) const {
    return (static_cast<size_t>(r) * static_cast<size_t>(cols)) + static_cast<size_t>(c);
  }
};

using InType = BinaryImage;
using OutType = HullList;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kamalagin_a_binary_image_convex_hull_omp
