#include "tsarkov_k_jarvis_convex_hull/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <thread>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "tsarkov_k_jarvis_convex_hull/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tsarkov_k_jarvis_convex_hull {

namespace {

struct BestPointState {
  std::size_t point_index = 0;
  bool is_initialized = false;
};

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

BestPointState FindBestStateInRange(const std::vector<Point> &unique_points, std::size_t current_point_index,
                                    std::size_t range_begin, std::size_t range_end) {
  BestPointState local_best_state;

  for (std::size_t point_index = range_begin; point_index < range_end; ++point_index) {
    if (point_index == current_point_index) {
      continue;
    }

    const BestPointState candidate_state{.point_index = point_index, .is_initialized = true};
    local_best_state = SelectBetterState(unique_points, current_point_index, local_best_state, candidate_state);
  }

  return local_best_state;
}

BestPointState FindBestStateOMP(const std::vector<Point> &unique_points, std::size_t current_point_index,
                                std::size_t range_begin, std::size_t range_end) {
  const auto range_size = range_end - range_begin;

  if (range_size == 0) {
    return BestPointState{};
  }

  const auto requested_thread_count = static_cast<std::size_t>(ppc::util::GetNumThreads());
  const auto actual_thread_count = std::min(requested_thread_count, range_size);
  const auto actual_thread_count_int = static_cast<int>(actual_thread_count);
  const auto chunk_size = (range_size + actual_thread_count - 1) / actual_thread_count;

  std::vector<BestPointState> local_best_states(actual_thread_count);

#pragma omp parallel for default(none) schedule(static) num_threads(actual_thread_count_int)          \
    shared(unique_points, current_point_index, range_begin, range_end, chunk_size, local_best_states, \
               actual_thread_count_int)
  for (int thread_index = 0; thread_index < actual_thread_count_int; ++thread_index) {
    const auto current_thread_index = static_cast<std::size_t>(thread_index);
    const std::size_t thread_range_begin = range_begin + (current_thread_index * chunk_size);
    const std::size_t thread_range_end = std::min(thread_range_begin + chunk_size, range_end);

    local_best_states[current_thread_index] =
        FindBestStateInRange(unique_points, current_point_index, thread_range_begin, thread_range_end);
  }

  BestPointState best_state;

  for (const BestPointState &local_best_state : local_best_states) {
    best_state = SelectBetterState(unique_points, current_point_index, best_state, local_best_state);
  }

  return best_state;
}

BestPointState FindBestStateSTL(const std::vector<Point> &unique_points, std::size_t current_point_index,
                                std::size_t range_begin, std::size_t range_end) {
  const auto range_size = range_end - range_begin;

  if (range_size == 0) {
    return BestPointState{};
  }

  const auto requested_thread_count = static_cast<std::size_t>(ppc::util::GetNumThreads());
  const auto actual_thread_count = std::min(requested_thread_count, range_size);
  const auto chunk_size = (range_size + actual_thread_count - 1) / actual_thread_count;

  std::vector<std::thread> worker_threads(actual_thread_count);
  std::vector<BestPointState> local_best_states(actual_thread_count);

  for (std::size_t thread_index = 0; thread_index < actual_thread_count; ++thread_index) {
    const std::size_t thread_range_begin = range_begin + (thread_index * chunk_size);
    const std::size_t thread_range_end = std::min(thread_range_begin + chunk_size, range_end);

    worker_threads[thread_index] = std::thread([&unique_points, &local_best_states, current_point_index,
                                                thread_range_begin, thread_range_end, thread_index]() {
      local_best_states[thread_index] =
          FindBestStateInRange(unique_points, current_point_index, thread_range_begin, thread_range_end);
    });
  }

  for (std::thread &worker_thread : worker_threads) {
    if (worker_thread.joinable()) {
      worker_thread.join();
    }
  }

  BestPointState best_state;

  for (const BestPointState &local_best_state : local_best_states) {
    best_state = SelectBetterState(unique_points, current_point_index, best_state, local_best_state);
  }

  return best_state;
}

BestPointState FindBestStateTBB(const std::vector<Point> &unique_points, std::size_t current_point_index,
                                std::size_t range_begin, std::size_t range_end) {
  if (range_begin == range_end) {
    return BestPointState{};
  }

  return oneapi::tbb::parallel_reduce(
      oneapi::tbb::blocked_range<std::size_t>(range_begin, range_end), BestPointState{},
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

std::size_t FindNextHullPointIndexALL(const std::vector<Point> &unique_points, std::size_t current_point_index) {
  const auto point_count = unique_points.size();
  const auto first_border = point_count / 3;
  const auto second_border = (2 * point_count) / 3;

  const BestPointState omp_best_state = FindBestStateOMP(unique_points, current_point_index, 0, first_border);
  const BestPointState stl_best_state =
      FindBestStateSTL(unique_points, current_point_index, first_border, second_border);
  const BestPointState tbb_best_state =
      FindBestStateTBB(unique_points, current_point_index, second_border, point_count);

  BestPointState result_state;
  result_state = SelectBetterState(unique_points, current_point_index, result_state, omp_best_state);
  result_state = SelectBetterState(unique_points, current_point_index, result_state, stl_best_state);
  result_state = SelectBetterState(unique_points, current_point_index, result_state, tbb_best_state);

  return result_state.point_index;
}

}  // namespace

TsarkovKJarvisConvexHullALL::TsarkovKJarvisConvexHullALL(const InType &input_points) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_points;
  GetOutput().clear();
}

bool TsarkovKJarvisConvexHullALL::ValidationImpl() {
  return !GetInput().empty() && GetOutput().empty();
}

bool TsarkovKJarvisConvexHullALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool TsarkovKJarvisConvexHullALL::RunImpl() {
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::vector<Point> unique_points = RemoveDuplicatePoints(GetInput());

  if (unique_points.empty()) {
    return false;
  }

  if (unique_points.size() == 1) {
    GetOutput() = unique_points;
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
  }

  if (unique_points.size() == 2) {
    GetOutput() = unique_points;
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
  }

  const auto start_point_index = FindLeftmostPointIndex(unique_points);
  std::size_t current_point_index = start_point_index;

  while (true) {
    GetOutput().push_back(unique_points[current_point_index]);

    current_point_index = FindNextHullPointIndexALL(unique_points, current_point_index);

    if (current_point_index == start_point_index) {
      break;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  return (rank >= 0) && !GetOutput().empty();
}

bool TsarkovKJarvisConvexHullALL::PostProcessingImpl() {
  return true;
}

}  // namespace tsarkov_k_jarvis_convex_hull
