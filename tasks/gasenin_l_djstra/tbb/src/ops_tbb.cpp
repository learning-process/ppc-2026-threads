#include "gasenin_l_djstra/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"

namespace gasenin_l_djstra {

GaseninLDjstraTBB::GaseninLDjstraTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool GaseninLDjstraTBB::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool GaseninLDjstraTBB::PreProcessingImpl() {
  vertex_count_ = GetInput();
  adj_.resize(vertex_count_);
  for (int i = 0; i < vertex_count_ - 1; ++i) {
    adj_[i].emplace_back(i + 1, 1);
  }
  dist_.assign(vertex_count_, kInf);
  dist_[0] = 0;
  return true;
}

bool GaseninLDjstraTBB::RunImpl() {
  std::vector<bool> visited(vertex_count_, false);

  for (int i = 0; i < vertex_count_; ++i) {
    int u = -1;
    int min_dist = kInf;
    for (int j = 0; j < vertex_count_; ++j) {
      if (!visited[j] && dist_[j] < min_dist) {
        min_dist = dist_[j];
        u = j;
      }
    }
    if (u == -1 || dist_[u] == kInf) {
      break;
    }
    visited[u] = true;

    const auto &neighbors = adj_[u];
    tbb::parallel_for(static_cast<size_t>(0), neighbors.size(), [&](size_t idx) {
      int v = neighbors[idx].first;
      int w = neighbors[idx].second;
      int new_dist = dist_[u] + w;
      dist_[v] = std::min(new_dist, dist_[v]);
    });
  }
  return true;
}

bool GaseninLDjstraTBB::PostProcessingImpl() {
  int sum = 0;
  for (int d : dist_) {
    if (d < kInf) {
      sum += d;
    }
  }
  GetOutput() = sum;
  return true;
}

}  // namespace gasenin_l_djstra
