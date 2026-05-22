#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "util/include/util.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

using VertexChoice = std::pair<Cost, NodeId>;

VertexChoice SelectBestVertex(const std::vector<char> &visited, const std::vector<Cost> &dist, int n, int threads) {
  Cost best_cost = kInf;
  NodeId best_vtx = -1;

#pragma omp parallel default(none) num_threads(threads) shared(visited, dist, n, best_cost, best_vtx)
  {
    Cost thr_cost = kInf;
    NodeId thr_vtx = -1;

#pragma omp for nowait
    for (int vi = 0; vi < n; ++vi) {
      const auto idx = static_cast<std::size_t>(vi);
      if (visited[idx] == 0 && dist[idx] < thr_cost) {
        thr_cost = dist[idx];
        thr_vtx = vi;
      }
    }

#pragma omp critical
    {
      if (thr_cost < best_cost || (thr_cost == best_cost && thr_vtx != -1 && (best_vtx == -1 || thr_vtx < best_vtx))) {
        best_cost = thr_cost;
        best_vtx = thr_vtx;
      }
    }
  }

  return {best_cost, best_vtx};
}

void RelaxNeighbors(const std::vector<std::pair<NodeId, Cost>> &nbrs, Cost best_cost, const std::vector<char> &visited,
                    std::vector<Cost> &dist, int threads) {
  const std::size_t sz = nbrs.size();
#pragma omp parallel for default(none) num_threads(threads) shared(nbrs, best_cost, visited, dist, sz) schedule(static)
  for (std::size_t i = 0; i < sz; ++i) {
    const NodeId tgt = nbrs[i].first;
    const Cost w = nbrs[i].second;
    const auto tgt_idx = static_cast<std::size_t>(tgt);

    if (visited[tgt_idx] == 0 && best_cost <= kInf - w) {
      const Cost cand = best_cost + w;
      std::atomic_ref<Cost> target(dist[tgt_idx]);
      Cost old_val = target.load(std::memory_order_relaxed);

      while (cand < old_val) {
        if (target.compare_exchange_weak(old_val, cand, std::memory_order_relaxed)) {
          break;
        }
      }
    }
  }
}

std::int64_t SumDistances(const std::vector<Cost> &dist, int n, int threads) {
  std::int64_t total = 0;

#pragma omp parallel for default(none) num_threads(threads) shared(dist, n) reduction(+ : total)
  for (int vi = 0; vi < n; ++vi) {
    const auto d = dist[static_cast<std::size_t>(vi)];
    if (d != kInf) {
      total += static_cast<std::int64_t>(d);
    }
  }

  return total;
}

}  // namespace

NalitovDDijkstrasAlgorithmOmp::NalitovDDijkstrasAlgorithmOmp(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmOmp::ValidationImpl() {
  if (GetOutput() != 0) {
    return false;
  }
  const InType &in = GetInput();
  if (in.n <= 0 || in.n > 10000) {
    return false;
  }
  if (in.source < 0 || in.source >= in.n) {
    return false;
  }

  return std::ranges::all_of(in.arcs, [&in](const Arc &a) {
    return a.from >= 0 && a.to >= 0 && a.from < in.n && a.to < in.n && a.weight >= 0;
  });
}

bool NalitovDDijkstrasAlgorithmOmp::PreProcessingImpl() {
  const InType &in = GetInput();
  const int vn = in.n;
  if (vn <= 0) {
    return false;
  }

  omp_threads_ = ppc::util::GetNumThreads();
  if (omp_threads_ <= 0) {
    return false;
  }

  const auto vn_u = static_cast<std::size_t>(vn);
  std::vector<std::size_t> outdeg(vn_u, 0);
  for (const Arc &a : in.arcs) {
    ++outdeg[static_cast<std::size_t>(a.from)];
  }

  graph_.assign(vn_u, {});
  std::vector<std::atomic<std::size_t>> row_next(vn_u);
  for (std::size_t i = 0; i < vn_u; ++i) {
    graph_[i].resize(outdeg[i]);
    row_next[i].store(0, std::memory_order_relaxed);
  }

  const std::size_t arc_count = in.arcs.size();
  auto &graph = graph_;

#pragma omp parallel for default(none) num_threads(omp_threads_) shared(in, row_next, graph, arc_count) schedule(guided)
  for (std::size_t ei = 0; ei < arc_count; ++ei) {
    const Arc &a = in.arcs[ei];
    const auto u = static_cast<std::size_t>(a.from);
    const std::size_t slot = row_next[u].fetch_add(1, std::memory_order_relaxed);
    graph[u][slot] = {a.to, a.weight};
  }

  dist_.assign(vn_u, kInf);
  visited_.assign(vn_u, 0);
  dist_[static_cast<std::size_t>(in.source)] = 0;

  return true;
}

bool NalitovDDijkstrasAlgorithmOmp::RunImpl() {
  const int n = GetInput().n;
  if (graph_.size() != static_cast<std::size_t>(n)) {
    return false;
  }

  auto &dist = dist_;
  auto &visited = visited_;
  const auto &graph = graph_;
  const int threads = omp_threads_;

  for (int step = 0; step < n; ++step) {
    const auto choice = SelectBestVertex(visited, dist, n, threads);
    const Cost best_cost = choice.first;
    const NodeId best_vtx = choice.second;

    if (best_vtx == -1 || best_cost == kInf) {
      break;
    }

    visited[static_cast<std::size_t>(best_vtx)] = 1;
    RelaxNeighbors(graph[static_cast<std::size_t>(best_vtx)], best_cost, visited, dist, threads);
  }

  GetOutput() = static_cast<OutType>(SumDistances(dist, n, threads));
  return true;
}

bool NalitovDDijkstrasAlgorithmOmp::PostProcessingImpl() {
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm
