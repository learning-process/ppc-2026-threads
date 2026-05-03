#include "tsarkov_k_jarvis_convex_hull/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <vector>

#include "tsarkov_k_jarvis_convex_hull/common/include/common.hpp"

namespace tsarkov_k_jarvis_convex_hull {

namespace {

std::int64_t CrossProduct(const Point &first_point, const Point &second_point, const Point &third_point) {
  const std::int64_t vector_first_x =
      static_cast<std::int64_t>(second_point.x) - static_cast<std::int64_t>(first_point.x);
  const std::int64_t vector_first_y =
      static_cast<std::int64_t>(second_point.y) - static_cast<std::int64_t>(first_point.y);
  const std::int64_t vector_second_x =
      static_cast<std::int64_t>(third_point.x) - static_cast<std::int64_t>(first_point.x);
  const std::int64_t vector_second_y =
      static_cast<std::int64_t>(third_point.y) - static_cast<std::int64_t>(first_point.y);

  return (vector_first_x * vector_second_y) - (vector_first_y * vector_second_x);
}

std::int64_t SquaredDistance(const Point &first_point, const Point &second_point) {
  const std::int64_t delta_x = static_cast<std::int64_t>(second_point.x) - static_cast<std::int64_t>(first_point.x);
  const std::int64_t delta_y = static_cast<std::int64_t>(second_point.y) - static_cast<std::int64_t>(first_point.y);

  return (delta_x * delta_x) + (delta_y * delta_y);
}

bool PointLess(const Point &first_point, const Point &second_point) {
  if (first_point.x != second_point.x) {
    return first_point.x < second_point.x;
  }
  return first_point.y < second_point.y;
}

std::vector<Point> RemoveDuplicatePoints(const std::vector<Point> &input_points) {
  std::vector<Point> unique_points = input_points;

  std::ranges::sort(unique_points, PointLess);
  unique_points.erase(std::ranges::unique(unique_points).begin(), unique_points.end());

  return unique_points;
}

std::size_t FindLeftmostPointIndex(const std::vector<Point> &input_points) {
  std::size_t leftmost_point_index = 0;

  for (std::size_t point_index = 1; point_index < input_points.size(); ++point_index) {
    const Point &current_point = input_points[point_index];
    const Point &leftmost_point = input_points[leftmost_point_index];

    if ((current_point.x < leftmost_point.x) ||
        ((current_point.x == leftmost_point.x) && (current_point.y < leftmost_point.y))) {
      leftmost_point_index = point_index;
    }
  }

  return leftmost_point_index;
}

std::size_t FindNextHullPointIndex(const std::vector<Point> &unique_points, std::size_t current_point_index) {
  std::size_t next_point_index = (current_point_index == 0) ? 1 : 0;

  for (std::size_t point_index = 0; point_index < unique_points.size(); ++point_index) {
    if (point_index == current_point_index) {
      continue;
    }

    const std::int64_t orientation =
        CrossProduct(unique_points[current_point_index], unique_points[next_point_index], unique_points[point_index]);

    if (orientation < 0) {
      next_point_index = point_index;
    } else if (orientation == 0) {
      const std::int64_t current_distance =
          SquaredDistance(unique_points[current_point_index], unique_points[next_point_index]);
      const std::int64_t candidate_distance =
          SquaredDistance(unique_points[current_point_index], unique_points[point_index]);

      if (candidate_distance > current_distance) {
        next_point_index = point_index;
      }
    }
  }

  return next_point_index;
}

}  // namespace

TsarkovKJarvisConvexHullSEQ::TsarkovKJarvisConvexHullSEQ(const InType &input_points) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_points;
  GetOutput().clear();
}

bool TsarkovKJarvisConvexHullSEQ::ValidationImpl() {
  return !GetInput().empty() && GetOutput().empty();
}

bool TsarkovKJarvisConvexHullSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool TsarkovKJarvisConvexHullSEQ::RunImpl() {
  const std::vector<Point> unique_points = RemoveDuplicatePoints(GetInput());

  if (unique_points.empty()) {
    return false;
  }

  if (unique_points.size() == 1) {
    GetOutput() = unique_points;
    return true;
  }

  if (unique_points.size() == 2) {
    GetOutput() = unique_points;
    return true;
  }

  const std::size_t start_point_index = FindLeftmostPointIndex(unique_points);
  std::size_t current_point_index = start_point_index;

  while (true) {
    GetOutput().push_back(unique_points[current_point_index]);

    current_point_index = FindNextHullPointIndex(unique_points, current_point_index);

    if (current_point_index == start_point_index) {
      break;
    }
  }

  return !GetOutput().empty();
}

bool TsarkovKJarvisConvexHullSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace tsarkov_k_jarvis_convex_hull
