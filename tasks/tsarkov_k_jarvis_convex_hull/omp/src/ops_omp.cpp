#include "tsarkov_k_jarvis_convex_hull/omp/include/ops_omp.hpp"

#include <omp.h>

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

bool ShouldReplaceBestPoint(const std::vector<Point> &unique_points, std::size_t current_point_index,
                            std::size_t best_point_index, std::size_t candidate_point_index) {
  if (candidate_point_index == current_point_index) {
    return false;
  }

  const std::int64_t orientation = CrossProduct(unique_points[current_point_index], unique_points[best_point_index],
                                                unique_points[candidate_point_index]);

  if (orientation < 0) {
    return true;
  }

  if (orientation == 0) {
    const std::int64_t best_distance =
        SquaredDistance(unique_points[current_point_index], unique_points[best_point_index]);
    const std::int64_t candidate_distance =
        SquaredDistance(unique_points[current_point_index], unique_points[candidate_point_index]);

    return candidate_distance > best_distance;
  }

  return false;
}

std::size_t FindThreadLocalBestIndex(const std::vector<Point> &unique_points, std::size_t current_point_index,
                                     std::size_t initial_best_index, std::size_t range_begin, std::size_t range_end) {
  std::size_t local_best_index = initial_best_index;

  for (std::size_t point_index = range_begin; point_index < range_end; ++point_index) {
    if (ShouldReplaceBestPoint(unique_points, current_point_index, local_best_index, point_index)) {
      local_best_index = point_index;
    }
  }

  return local_best_index;
}

std::size_t ReduceBestIndices(const std::vector<Point> &unique_points, std::size_t current_point_index,
                              std::size_t initial_best_index, const std::vector<std::size_t> &local_best_indices) {
  std::size_t next_point_index = initial_best_index;

  for (std::size_t local_best_index : local_best_indices) {
    if (ShouldReplaceBestPoint(unique_points, current_point_index, next_point_index, local_best_index)) {
      next_point_index = local_best_index;
    }
  }

  return next_point_index;
}

std::size_t FindNextHullPointIndexOMP(const std::vector<Point> &unique_points, std::size_t current_point_index) {
  const std::size_t initial_next_point_index = (current_point_index == 0) ? 1 : 0;
  const int thread_count = omp_get_max_threads();

  std::vector<std::size_t> local_best_indices(static_cast<std::size_t>(thread_count), initial_next_point_index);

#pragma omp parallel default(none) \
    shared(unique_points, current_point_index, local_best_indices, initial_next_point_index, thread_count)
  {
    const int thread_index = omp_get_thread_num();
    const std::size_t point_count = unique_points.size();
    const std::size_t chunk_size =
        (point_count + static_cast<std::size_t>(thread_count) - 1) / static_cast<std::size_t>(thread_count);
    const std::size_t range_begin = static_cast<std::size_t>(thread_index) * chunk_size;
    const std::size_t range_end = std::min(range_begin + chunk_size, point_count);

    local_best_indices[static_cast<std::size_t>(thread_index)] =
        FindThreadLocalBestIndex(unique_points, current_point_index, initial_next_point_index, range_begin, range_end);
  }

  return ReduceBestIndices(unique_points, current_point_index, initial_next_point_index, local_best_indices);
}

}  // namespace

TsarkovKJarvisConvexHullOMP::TsarkovKJarvisConvexHullOMP(const InType &input_points) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_points;
  GetOutput().clear();
}

bool TsarkovKJarvisConvexHullOMP::ValidationImpl() {
  return !GetInput().empty() && GetOutput().empty();
}

bool TsarkovKJarvisConvexHullOMP::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool TsarkovKJarvisConvexHullOMP::RunImpl() {
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

    current_point_index = FindNextHullPointIndexOMP(unique_points, current_point_index);

    if (current_point_index == start_point_index) {
      break;
    }
  }

  return !GetOutput().empty();
}

bool TsarkovKJarvisConvexHullOMP::PostProcessingImpl() {
  return true;
}

}  // namespace tsarkov_k_jarvis_convex_hull
