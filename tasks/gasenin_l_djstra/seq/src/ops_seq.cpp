#include "gasenin_l_djstra/seq/include/ops_seq.hpp"

#include <limits>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"

namespace gasenin_l_djstra {

GaseninLDjstraSEQ::GaseninLDjstraSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool GaseninLDjstraSEQ::ValidationImpl() {
  return GetInput() > 0;
}

bool GaseninLDjstraSEQ::PreProcessingImpl() {
  return true;
}

bool GaseninLDjstraSEQ::RunImpl() {
  InType n = GetInput();
  const InType INF = std::numeric_limits<InType>::max();
  std::vector<InType> dist(n, INF);
  std::vector<bool> visited(n, false);
  dist[0] = 0;

  for (int count = 0; count < n; ++count) {
    InType u = -1;
    InType min_dist = INF;
    for (int i = 0; i < n; ++i) {
      if (!visited[i] && dist[i] < min_dist) {
        min_dist = dist[i];
        u = i;
      }
    }
    if (u == -1) {
      break;
    }
    visited[u] = true;

    for (int v = 0; v < n; ++v) {
      if (!visited[v] && u != v) {
        InType weight = (u > v) ? (u - v) : (v - u);
        if (dist[u] != INF && dist[u] + weight < dist[v]) {
          dist[v] = dist[u] + weight;
        }
      }
    }
  }

  InType sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += dist[i];
  }
  GetOutput() = sum;
  return true;
}

bool GaseninLDjstraSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace gasenin_l_djstra
