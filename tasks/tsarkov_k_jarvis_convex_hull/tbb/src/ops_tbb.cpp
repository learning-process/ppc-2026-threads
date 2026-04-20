#include "tsarkov_k_jarvis_convex_hull/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_reduce.h>

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

struct BestPointState {
  std::size_t point_index = 0;
  bool is_initialized = false;
};

BestPointState SelectBetterState(const std::vector<Point> &unique_points, std::size_t current_point_index,
                                 const BestPointState &first_state, const BestPointState &second_state) {
  if (!first_state.is_initialized) {
    return second_state;
  }
  if (!second_state.is_initialized) {
    return first_state;
  }

  if (ShouldReplaceBestPoint(unique_points, current_point_index, first_state.point_index, second_state.point_index)) {
    return second_state;
  }

  return first_state;
}

BestPointState FindNextHullPointStateTBB(const std::vector<Point> &unique_points, std::size_t current_point_index) {
  return oneapi::tbb::parallel_reduce(
      oneapi::tbb::blocked_range<std::size_t>(0, unique_points.size()), BestPointState{},
      [&](const oneapi::tbb::blocked_range<std::size_t> &point_range, BestPointState local_best_state) {
    for (std::size_t point_index = point_range.begin(); point_index < point_range.end(); ++point_index) {
      if (point_index == current_point_index) {
        continue;
      }

      const BestPointState candidate_state{.point_index = point_index, .is_initialized = true};
      local_best_state = SelectBetterState(unique_points, current_point_index, local_best_state, candidate_state);
    }

    return local_best_state;
  }, [&](const BestPointState &first_state, const BestPointState &second_state) {
    return SelectBetterState(unique_points, current_point_index, first_state, second_state);
  });
}

std::size_t FindNextHullPointIndexTBB(const std::vector<Point> &unique_points, std::size_t current_point_index) {
  const BestPointState best_point_state = FindNextHullPointStateTBB(unique_points, current_point_index);
  return best_point_state.point_index;
}

}  // namespace

TsarkovKJarvisConvexHullTBB::TsarkovKJarvisConvexHullTBB(const InType &input_points) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_points;
  GetOutput().clear();
}

bool TsarkovKJarvisConvexHullTBB::ValidationImpl() {
  return !GetInput().empty() && GetOutput().empty();
}

bool TsarkovKJarvisConvexHullTBB::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool TsarkovKJarvisConvexHullTBB::RunImpl() {
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

  const auto start_point_index = FindLeftmostPointIndex(unique_points);
  std::size_t current_point_index = start_point_index;

  while (true) {
    GetOutput().push_back(unique_points[current_point_index]);

    current_point_index = FindNextHullPointIndexTBB(unique_points, current_point_index);

    if (current_point_index == start_point_index) {
      break;
    }
  }

  return !GetOutput().empty();
}

bool TsarkovKJarvisConvexHullTBB::PostProcessingImpl() {
  return true;
}

}  // namespace tsarkov_k_jarvis_convex_hull
