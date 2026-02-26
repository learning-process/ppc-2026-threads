#include "egorova_l_binary_convex_hull/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stack>
#include <tuple>
#include <utility>  // для std::move
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"

namespace {

using Point = egorova_l_binary_convex_hull::Point;

int64_t CrossProduct(const Point &a, const Point &b, const Point &c) {
  return (static_cast<int64_t>(b.x - a.x) * (c.y - a.y)) - (static_cast<int64_t>(b.y - a.y) * (c.x - a.x));
}

// Helper function for building lower hull
void BuildLowerHull(const std::vector<Point> &points, size_t n, std::vector<Point> &hull) {
  for (size_t i = 0; i < n; ++i) {
    while (hull.size() >= 2) {
      const Point &a = hull[hull.size() - 2];
      const Point &b = hull.back();
      if (CrossProduct(a, b, points[i]) <= 0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(points[i]);
  }
}

// Helper function for building upper hull
void BuildUpperHull(const std::vector<Point> &points, size_t n, size_t lower_size, std::vector<Point> &hull) {
  for (size_t i = n - 1; i > 0; --i) {
    const size_t idx = i - 1;
    while (hull.size() > lower_size) {
      const Point &a = hull[hull.size() - 2];
      const Point &b = hull.back();
      if (CrossProduct(a, b, points[idx]) <= 0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(points[idx]);
  }
}

void BuildConvexHull(std::vector<Point> &points, std::vector<Point> &hull) {
  const size_t n = points.size();
  if (n <= 2) {
    hull.assign(points.begin(), points.end());
    return;
  }

  // Sort points in-place
  std::ranges::sort(points,
                    [](const Point &lhs, const Point &rhs) { return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y); });

  hull.clear();
  hull.reserve(n + 1);  // Reserve space to avoid reallocations

  // Build lower hull
  BuildLowerHull(points, n, hull);

  // Build upper hull
  const size_t lower_size = hull.size();
  BuildUpperHull(points, n, lower_size, hull);

  // Remove the last point if it's the same as the first (closed polygon)
  if (hull.size() > 1 && hull.front().x == hull.back().x && hull.front().y == hull.back().y) {
    hull.pop_back();
  }
}

// Helper function to process a single direction
inline void TryAddNeighbor(const std::vector<uint8_t> &image, int width, int height, int next_x, int next_y, int label,
                           std::vector<int> &labels, std::stack<Point> &stack) {
  if (next_x >= 0 && next_x < width && next_y >= 0 && next_y < height) {
    const size_t next_index = (static_cast<size_t>(next_y) * static_cast<size_t>(width)) + static_cast<size_t>(next_x);
    if (image[next_index] != 0 && labels[next_index] == 0) {
      labels[next_index] = label;
      stack.push({next_x, next_y});
    }
  }
}

void ProcessComponent(const std::vector<uint8_t> &image, int width, int height, int start_x, int start_y, int label,
                      std::vector<int> &labels, std::vector<Point> &component_points) {
  std::stack<Point> stack;
  stack.push({start_x, start_y});

  const size_t start_index = (static_cast<size_t>(start_y) * static_cast<size_t>(width)) + static_cast<size_t>(start_x);
  labels[start_index] = label;

  component_points.clear();
  component_points.reserve(100);  // Reserve some space to avoid reallocations

  while (!stack.empty()) {
    const Point current = stack.top();
    stack.pop();
    component_points.push_back(current);

    // Process all 4 directions using helper function
    TryAddNeighbor(image, width, height, current.x + 1, current.y, label, labels, stack);
    TryAddNeighbor(image, width, height, current.x - 1, current.y, label, labels, stack);
    TryAddNeighbor(image, width, height, current.x, current.y + 1, label, labels, stack);
    TryAddNeighbor(image, width, height, current.x, current.y - 1, label, labels, stack);
  }
}

}  // namespace

namespace egorova_l_binary_convex_hull {

BinaryConvexHullSEQ::BinaryConvexHullSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BinaryConvexHullSEQ::ValidationImpl() {
  const auto &input = GetInput();
  return input.width > 0 && input.height > 0 && !input.data.empty();
}

bool BinaryConvexHullSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool BinaryConvexHullSEQ::RunImpl() {
  const int width = GetInput().width;
  const int height = GetInput().height;
  const auto &image = GetInput().data;
  const size_t image_size = static_cast<size_t>(width) * static_cast<size_t>(height);

  std::vector<int> labels(image_size, 0);
  int label_counter = 0;

  // Clear output and reserve some space
  auto &output = GetOutput();
  output.clear();
  output.reserve(10);  // Reserve space for expected number of components

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      const size_t index = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);
      if (image[index] != 0 && labels[index] == 0) {
        ++label_counter;
        std::vector<Point> component_points;

        ProcessComponent(image, width, height, col, row, label_counter, labels, component_points);

        if (!component_points.empty()) {
          std::vector<Point> hull;
          BuildConvexHull(component_points, hull);
          output.push_back(std::move(hull));
        }
      }
    }
  }
  return true;
}

bool BinaryConvexHullSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace egorova_l_binary_convex_hull
