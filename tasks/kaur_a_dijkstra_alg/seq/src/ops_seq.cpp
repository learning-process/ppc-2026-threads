#include "../include/seq.hpp"

#include <limits>
#include <vector>

#include "../../common/include/common.hpp"

namespace kaur_a_dijkstra_alg {

KaurADijkstraAlgSEQ::KaurADijkstraAlgSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KaurADijkstraAlgSEQ::ValidationImpl() {
  return GetInput() > 0;
}

bool KaurADijkstraAlgSEQ::PreProcessingImpl() {
  return true;
}

// static
InType KaurADijkstraAlgSEQ::FindMinDist(const std::vector<InType> &dist, const std::vector<bool> &visited) {
  auto n = static_cast<InType>(dist.size());
  InType min_dist = std::numeric_limits<InType>::max();
  InType u = -1;
  for (InType i = 0; i < n; ++i) {
    if (!visited[i] && dist[i] < min_dist) {
      min_dist = dist[i];
      u = i;
    }
  }
  return u;
}

// static
void KaurADijkstraAlgSEQ::RelaxEdges(InType u, std::vector<InType> &dist, const std::vector<bool> &visited) {
  auto n = static_cast<InType>(dist.size());
  for (InType vert = 0; vert < n; ++vert) {
    if (!visited[vert] && u != vert) {
      InType weight = (u > vert) ? (u - vert) : (vert - u);
      if (dist[u] != std::numeric_limits<InType>::max() && dist[u] + weight < dist[vert]) {
        dist[vert] = dist[u] + weight;
      }
    }
  }
}

bool KaurADijkstraAlgSEQ::RunImpl() {
  InType n = GetInput();
  if (n == 0) {
    return false;
  }

  const InType k_inf = std::numeric_limits<InType>::max();
  std::vector<InType> dist(n, k_inf);
  std::vector<bool> visited(n, false);
  dist[0] = 0;

  for (int count = 0; count < n; ++count) {
    InType u = FindMinDist(dist, visited);
    if (u == -1) {
      break;
    }
    visited[u] = true;
    RelaxEdges(u, dist, visited);
  }

  InType sum = 0;
  for (InType i = 0; i < n; ++i) {
    sum += dist[i];
  }
  GetOutput() = sum;
  return true;
}

bool KaurADijkstraAlgSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kaur_a_dijkstra_alg
