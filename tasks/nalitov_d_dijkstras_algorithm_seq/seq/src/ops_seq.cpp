#include "nalitov_d_dijkstras_algorithm_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "nalitov_d_dijkstras_algorithm_seq/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm_seq {

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

InType FindMinDistanceVertex(const std::vector<InType> &distances, const std::vector<bool> &processed,
                             InType k_infinity) {
  InType min_distance = k_infinity;
  InType current_vertex = -1;
  const auto num_vertices = static_cast<InType>(distances.size());

  for (InType vertex_idx = 0; vertex_idx < num_vertices; ++vertex_idx) {
    if (!processed[vertex_idx] && distances[vertex_idx] < min_distance) {
      min_distance = distances[vertex_idx];
      current_vertex = vertex_idx;
    }
  }

  return current_vertex;
}

void RelaxEdges(InType current_vertex, std::vector<InType> &distances, const std::vector<bool> &processed,
                InType k_infinity) {
  const auto num_vertices = static_cast<InType>(distances.size());

  for (InType neighbor = 0; neighbor < num_vertices; ++neighbor) {
    if (processed[neighbor] || neighbor == current_vertex) {
      continue;
    }

    const InType edge_weight = GetEdgeWeight(current_vertex, neighbor);

    if (distances[current_vertex] == k_infinity) {
      continue;
    }

    InType new_distance = 0;
    if (!SafeAdd(distances[current_vertex], edge_weight, new_distance)) {
      continue;
    }

    distances[neighbor] = std::min(new_distance, distances[neighbor]);
  }
}

OutType CalculateTotalDistance(const std::vector<InType> &distances, InType k_infinity) {
  OutType total_sum = 0;
  const auto num_vertices = static_cast<InType>(distances.size());

  for (InType vertex_idx = 0; vertex_idx < num_vertices; ++vertex_idx) {
    if (distances[vertex_idx] != k_infinity) {
      if (total_sum > std::numeric_limits<OutType>::max() - distances[vertex_idx]) {
        return -1;  // Overflow indicator
      }
      total_sum += distances[vertex_idx];
    }
  }

  return total_sum;
}

}  // namespace

NalitovDDijkstrasAlgorithmSeq::NalitovDDijkstrasAlgorithmSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmSeq::ValidationImpl() {
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

bool NalitovDDijkstrasAlgorithmSeq::PreProcessingImpl() {
  GetOutput() = 0;
  return true;
}

bool NalitovDDijkstrasAlgorithmSeq::RunImpl() {
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
  std::vector<InType> distances(n, k_infinity);
  std::vector<bool> processed(n, false);

  if (distances.empty()) {
    return false;
  }
  distances[0] = 0;

  for (InType iteration = 0; iteration < n; ++iteration) {
    const InType current_vertex = FindMinDistanceVertex(distances, processed, k_infinity);

    if (current_vertex == -1 || distances[current_vertex] == k_infinity) {
      break;
    }

    processed[current_vertex] = true;
    RelaxEdges(current_vertex, distances, processed, k_infinity);
  }

  const OutType total_sum = CalculateTotalDistance(distances, k_infinity);
  if (total_sum < 0) {
    return false;  // Overflow occurred
  }

  GetOutput() = total_sum;
  return true;
}

bool NalitovDDijkstrasAlgorithmSeq::PostProcessingImpl() {
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm_seq
