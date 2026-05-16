#include "nalitov_d_dijkstras_algorithm/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_reduce.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

struct BestVertex {
  Cost cost{kInf};
  NodeId vtx{-1};

  [[nodiscard]] bool IsBetterThan(const BestVertex &other) const {
    return cost < other.cost || (cost == other.cost && vtx != -1 && (other.vtx == -1 || vtx < other.vtx));
  }

  void Merge(const BestVertex &other) {
    if (other.IsBetterThan(*this)) {
      cost = other.cost;
      vtx = other.vtx;
    }
  }
};

BestVertex FindLocalBestVertex(const std::vector<char> &visited, const std::vector<std::atomic<Cost>> &atomic_dist,
                               const tbb::blocked_range<int> &range, BestVertex acc) {
  for (int vi = range.begin(); vi < range.end(); ++vi) {
    const auto idx = static_cast<std::size_t>(vi);

    if (visited[idx] != 0) {
      continue;
    }

    const Cost dist = atomic_dist[idx].load(std::memory_order_relaxed);

    const BestVertex candidate{
        .cost = dist,
        .vtx = vi,
    };

    if (candidate.IsBetterThan(acc)) {
      acc = candidate;
    }
  }

  return acc;
}

BestVertex SelectBestVertex(const std::vector<char> &visited, const std::vector<std::atomic<Cost>> &atomic_dist,
                            int n) {
  return tbb::parallel_reduce(tbb::blocked_range<int>(0, n), BestVertex{}, [&](const auto &range, BestVertex acc) {
    return FindLocalBestVertex(visited, atomic_dist, range, acc);
  }, [](BestVertex lhs, const BestVertex &rhs) {
    lhs.Merge(rhs);
    return lhs;
  });
}

void RelaxEdge(std::atomic<Cost> &target, Cost candidate) {
  Cost old = target.load(std::memory_order_relaxed);

  while (candidate < old) {
    if (target.compare_exchange_weak(old, candidate, std::memory_order_relaxed)) {
      break;
    }
  }
}

void RelaxNeighbors(const std::vector<std::pair<NodeId, Cost>> &neighbors, Cost pivot_cost,
                    const std::vector<char> &visited, std::vector<std::atomic<Cost>> &atomic_dist) {
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, neighbors.size()), [&](const auto &range) {
    for (std::size_t ei = range.begin(); ei < range.end(); ++ei) {
      const auto [target, weight] = neighbors[ei];
      const auto target_idx = static_cast<std::size_t>(target);

      if (visited[target_idx] != 0 || pivot_cost > kInf - weight) {
        continue;
      }

      RelaxEdge(atomic_dist[target_idx], pivot_cost + weight);
    }
  });
}

std::int64_t SumDistances(const std::vector<std::atomic<Cost>> &atomic_dist, std::size_t n) {
  return tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, n), std::int64_t{0},
                              [&](const auto &range, std::int64_t acc) {
    for (std::size_t i = range.begin(); i < range.end(); ++i) {
      const Cost dist = atomic_dist[i].load(std::memory_order_relaxed);

      if (dist != kInf) {
        acc += static_cast<std::int64_t>(dist);
      }
    }

    return acc;
  }, std::plus<>());
}

}  // namespace

NalitovDDijkstrasAlgorithmTBB::NalitovDDijkstrasAlgorithmTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmTBB::ValidationImpl() {
  const InType &in = GetInput();

  if (in.n <= 0 || in.n > 10000 || in.source < 0 || in.source >= in.n) {
    return false;
  }

  return std::ranges::all_of(in.arcs, [&in](const Arc &arc) {
    return arc.from >= 0 && arc.to >= 0 && arc.from < in.n && arc.to < in.n && arc.weight >= 0;
  });
}

bool NalitovDDijkstrasAlgorithmTBB::PreProcessingImpl() {
  const InType &in = GetInput();
  const auto vertex_count = static_cast<std::size_t>(in.n);

  std::vector<std::size_t> out_degree(vertex_count, 0);

  for (const Arc &arc : in.arcs) {
    ++out_degree[static_cast<std::size_t>(arc.from)];
  }

  graph_.assign(vertex_count, {});

  std::vector<std::atomic<std::size_t>> row_next(vertex_count);

  for (std::size_t vi = 0; vi < vertex_count; ++vi) {
    graph_[vi].resize(out_degree[vi]);
    row_next[vi].store(0, std::memory_order_relaxed);
  }

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, in.arcs.size()), [&](const auto &range) {
    for (std::size_t ei = range.begin(); ei < range.end(); ++ei) {
      const Arc &arc = in.arcs[ei];
      const auto from = static_cast<std::size_t>(arc.from);

      const std::size_t slot = row_next[from].fetch_add(1, std::memory_order_relaxed);
      graph_[from][slot] = {arc.to, arc.weight};
    }
  });

  dist_.assign(vertex_count, kInf);
  visited_.assign(vertex_count, 0);
  dist_[static_cast<std::size_t>(in.source)] = 0;

  return true;
}

bool NalitovDDijkstrasAlgorithmTBB::RunImpl() {
  const int n = GetInput().n;
  const auto vertex_count = static_cast<std::size_t>(n);

  std::vector<std::atomic<Cost>> atomic_dist(vertex_count);

  for (std::size_t vi = 0; vi < vertex_count; ++vi) {
    atomic_dist[vi].store(dist_[vi], std::memory_order_relaxed);
  }

  for (int step = 0; step < n; ++step) {
    const BestVertex pivot = SelectBestVertex(visited_, atomic_dist, n);

    if (pivot.vtx == -1 || pivot.cost == kInf) {
      break;
    }

    visited_[static_cast<std::size_t>(pivot.vtx)] = 1;

    const auto &neighbors = graph_[static_cast<std::size_t>(pivot.vtx)];
    RelaxNeighbors(neighbors, pivot.cost, visited_, atomic_dist);
  }

  GetOutput() = static_cast<OutType>(SumDistances(atomic_dist, vertex_count));

  return true;
}

bool NalitovDDijkstrasAlgorithmTBB::PostProcessingImpl() {
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm
