#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

inline bool SafeAdd(InType a, InType b, InType &result) {
  if (a > 0 && b > std::numeric_limits<InType>::max() - a) {
    return false;
  }
  if (a < 0 && b < std::numeric_limits<InType>::min() - a) {
    return false;
  }
  result = a + b;
  return true;
}

inline InType GetEdgeWeight(InType from, InType to) {
  if (from == to) {
    return 0;
  }
  return (from > to) ? (from - to) : (to - from);
}

int GetWorkerCount() {
  int worker_count = 1;
#pragma omp parallel default(none) shared(worker_count)
  {
#pragma omp single
    worker_count = omp_get_num_threads();
  }
  return worker_count;
}

InType SelectNextVertexOmp(std::vector<InType> &distances, std::vector<char> &processed, InType k_infinity,
                           int worker_count) {
  const auto vertex_count = static_cast<InType>(distances.size());
  std::vector<InType> thread_best_distance(worker_count, k_infinity);
  std::vector<InType> thread_best_vertex(worker_count, -1);
  InType selected_vertex = -1;

#pragma omp parallel default(none) shared(vertex_count, distances, processed, thread_best_distance, \
                                              thread_best_vertex, worker_count, selected_vertex, k_infinity)
  {
    const int tid = omp_get_thread_num();
    InType best_distance = std::numeric_limits<InType>::max();
    InType best_vertex = -1;

#pragma omp for nowait
    for (InType vertex_idx = 0; vertex_idx < vertex_count; ++vertex_idx) {
      if (processed[vertex_idx] == 0 && distances[vertex_idx] < best_distance) {
        best_distance = distances[vertex_idx];
        best_vertex = vertex_idx;
      }
    }

    thread_best_distance[tid] = best_distance;
    thread_best_vertex[tid] = best_vertex;

#pragma omp barrier
#pragma omp single
    {
      InType best_global_distance = std::numeric_limits<InType>::max();
      for (int thread_idx = 0; thread_idx < worker_count; ++thread_idx) {
        if (thread_best_vertex[thread_idx] != -1 && (thread_best_distance[thread_idx] < best_global_distance ||
                                                     (thread_best_distance[thread_idx] == best_global_distance &&
                                                      thread_best_vertex[thread_idx] < selected_vertex))) {
          best_global_distance = thread_best_distance[thread_idx];
          selected_vertex = thread_best_vertex[thread_idx];
        }
      }
      if (selected_vertex != -1 && best_global_distance != k_infinity) {
        processed[selected_vertex] = 1;
      } else {
        selected_vertex = -1;
      }
    }
  }
  return selected_vertex;
}

void UpdateNeighborsOmp(InType anchor_vertex, std::vector<InType> &distances, const std::vector<char> &processed,
                        InType k_infinity) {
  const auto vertex_count = static_cast<InType>(distances.size());
  const InType anchor_distance = distances[anchor_vertex];

#pragma omp parallel for default(none) \
    shared(anchor_vertex, distances, processed, vertex_count, k_infinity, anchor_distance)
  for (InType neighbor = 0; neighbor < vertex_count; ++neighbor) {
    if (processed[neighbor] == 0 && neighbor != anchor_vertex) {
      const InType edge_weight = GetEdgeWeight(anchor_vertex, neighbor);
      if (anchor_distance == k_infinity) {
        continue;
      }
      InType new_distance = 0;
      if (!SafeAdd(anchor_distance, edge_weight, new_distance)) {
        continue;
      }
      distances[neighbor] = std::min(new_distance, distances[neighbor]);
    }
  }
}

OutType AccumulateReachableDistanceOmp(const std::vector<InType> &distances, InType k_infinity) {
  int64_t aggregated_distance = 0;
  const auto vertex_count = static_cast<InType>(distances.size());
#pragma omp parallel for reduction(+ : aggregated_distance) default(none) shared(distances, vertex_count, k_infinity)
  for (InType vertex_idx = 0; vertex_idx < vertex_count; ++vertex_idx) {
    if (distances[vertex_idx] != k_infinity) {
      aggregated_distance += distances[vertex_idx];
    }
  }
  if (aggregated_distance < 0 || aggregated_distance > std::numeric_limits<OutType>::max()) {
    return -1;
  }
  return static_cast<OutType>(aggregated_distance);
}

}  // namespace

NalitovDDijkstrasAlgorithmOmp::NalitovDDijkstrasAlgorithmOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmOmp::ValidationImpl() {
  const InType n = GetInput();

  constexpr InType kMaxVertices = 10000;
  if (n <= 0 || n > kMaxVertices) {
    return false;
  }

  if (GetOutput() != 0) {
    return false;
  }

  return true;
}

bool NalitovDDijkstrasAlgorithmOmp::PreProcessingImpl() {
  const InType n = GetInput();
  const InType k_infinity = std::numeric_limits<InType>::max();
  distances_.assign(n, k_infinity);
  processed_.assign(n, 0);
  if (!distances_.empty()) {
    distances_[0] = 0;
  }
  GetOutput() = 0;
  return true;
}

bool NalitovDDijkstrasAlgorithmOmp::RunImpl() {
  const InType n = GetInput();

  if (n <= 0) {
    return false;
  }

  if (n == 1) {
    GetOutput() = 0;
    return true;
  }

  if (n < 2) {
    return false;
  }

  const InType k_infinity = std::numeric_limits<InType>::max();
  if (distances_.empty()) {
    return false;
  }
  const int worker_count = GetWorkerCount();

  for (InType iteration = 0; iteration < n; ++iteration) {
    const InType current_vertex = SelectNextVertexOmp(distances_, processed_, k_infinity, worker_count);

    if (current_vertex == -1 || distances_[current_vertex] == k_infinity) {
      break;
    }

    UpdateNeighborsOmp(current_vertex, distances_, processed_, k_infinity);
  }

  const OutType total_sum = AccumulateReachableDistanceOmp(distances_, k_infinity);
  if (total_sum < 0) {
    return false;
  }

  GetOutput() = total_sum;
  return true;
}

bool NalitovDDijkstrasAlgorithmOmp::PostProcessingImpl() {
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm
