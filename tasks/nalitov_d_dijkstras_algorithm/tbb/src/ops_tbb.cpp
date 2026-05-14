#include "nalitov_d_dijkstras_algorithm/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <mutex>
#include <queue>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;

struct QueueItem {
  Cost dist{};
  NodeId node{};
  friend bool operator<(const QueueItem &a, const QueueItem &b) {
    if (a.dist != b.dist) {
      return a.dist < b.dist;
    }
    return a.node < b.node;
  }
  friend bool operator>(const QueueItem &a, const QueueItem &b) {
    return b < a;
  }
};

using DijkstraHeap = std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<>>;

bool SafeAdd(std::int64_t acc, Cost value, std::int64_t &result) {
  const auto v = static_cast<std::int64_t>(value);
  if (v > 0 && acc > std::numeric_limits<std::int64_t>::max() - v) {
    return false;
  }
  if (v < 0 && acc < std::numeric_limits<std::int64_t>::min() - v) {
    return false;
  }
  result = acc + v;
  return true;
}

void BuildGraphParallel(const InType &input, OutgoingTable &graph) {
  const int vn = input.n;
  const std::size_t arc_count = input.arcs.size();

  std::vector<std::size_t> outdegree(static_cast<std::size_t>(vn), 0);
  for (const Arc &a : input.arcs) {
    ++outdegree[static_cast<std::size_t>(a.from)];
  }

  graph.assign(static_cast<std::size_t>(vn), {});
  for (int vx = 0; vx < vn; ++vx) {
    graph[static_cast<std::size_t>(vx)].resize(outdegree[static_cast<std::size_t>(vx)]);
  }

  const auto vn_u = static_cast<std::size_t>(vn);
  std::vector<std::atomic<std::size_t>> slot_counter(vn_u);
  for (std::size_t vx = 0; vx < vn_u; ++vx) {
    slot_counter[vx].store(0, std::memory_order_relaxed);
  }

  OutgoingTable &g = graph;
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, arc_count), [&](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t i = range.begin(); i < range.end(); ++i) {
      const Arc &arc = input.arcs[i];
      const auto from = static_cast<std::size_t>(arc.from);
      const std::size_t pos = slot_counter[from].fetch_add(1, std::memory_order_relaxed);
      g[from][pos] = {arc.to, arc.weight};
    }
  });
}

void RelaxEdgesParallel(Cost d_u, const std::vector<std::pair<NodeId, Cost>> &out_edges, std::vector<Cost> &distance,
                        std::vector<std::pair<NodeId, Cost>> &updates) {
  const auto edge_cnt = out_edges.size();
  std::mutex upd_mutex;

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edge_cnt), [&](const tbb::blocked_range<std::size_t> &r) {
    for (std::size_t j = r.begin(); j < r.end(); ++j) {
      const auto &[v, w] = out_edges[j];
      if (d_u <= kInf - w) {
        const Cost cand = d_u + w;
        if (cand < distance[v]) {
          std::scoped_lock lock(upd_mutex);
          updates.emplace_back(v, cand);
        }
      }
    }
  });
}

void ApplyUpdates(const std::vector<std::pair<NodeId, Cost>> &updates, std::vector<Cost> &distance,
                  DijkstraHeap &heap) {
  for (const auto &[v, cand] : updates) {
    if (cand < distance[v]) {
      distance[v] = cand;
      heap.push({cand, v});
    }
  }
}

std::vector<Cost> ComputeShortestPaths(NodeId source, const OutgoingTable &graph) {
  const auto vertex_count = static_cast<NodeId>(graph.size());
  std::vector<Cost> distance(vertex_count, kInf);
  DijkstraHeap heap;
  std::vector<char> visited(vertex_count, 0);

  distance[source] = 0;
  heap.push({0, source});

  while (!heap.empty()) {
    const auto [d_u, u] = heap.top();
    heap.pop();

    if (visited[u] != 0) {
      continue;
    }
    visited[u] = 1;

    const auto &out_edges = graph[u];
    if (out_edges.empty()) {
      continue;
    }

    std::vector<std::pair<NodeId, Cost>> updates;
    updates.reserve(out_edges.size());
    RelaxEdgesParallel(d_u, out_edges, distance, updates);
    ApplyUpdates(updates, distance, heap);
  }

  return distance;
}

bool SumReachableDistances(const std::vector<Cost> &distance, OutType &total) {
  std::int64_t accumulator = 0;
  for (Cost d : distance) {
    if (d == kInf) {
      continue;
    }
    if (!SafeAdd(accumulator, d, accumulator)) {
      return false;
    }
  }
  if (accumulator < 0 || accumulator > std::numeric_limits<OutType>::max()) {
    return false;
  }
  total = accumulator;
  return true;
}

}  // namespace

NalitovDDijkstrasAlgorithmTBB::NalitovDDijkstrasAlgorithmTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmTBB::ValidationImpl() {
  if (GetOutput() != 0) {
    return false;
  }

  const InType &in = GetInput();
  constexpr int kMaxVertices = 10000;
  if (in.n <= 0 || in.n > kMaxVertices) {
    return false;
  }
  if (in.source < 0 || in.source >= in.n) {
    return false;
  }

  const auto arc_valid = [&in](const Arc &a) {
    return a.from >= 0 && a.to >= 0 && a.from < in.n && a.to < in.n && a.weight >= 0;
  };
  return std::ranges::all_of(in.arcs, arc_valid);
}

bool NalitovDDijkstrasAlgorithmTBB::PreProcessingImpl() {
  const InType &in = GetInput();
  BuildGraphParallel(in, graph_);
  GetOutput() = 0;
  return true;
}

bool NalitovDDijkstrasAlgorithmTBB::RunImpl() {
  const InType &in = GetInput();
  if (graph_.size() != static_cast<std::size_t>(in.n)) {
    return false;
  }

  const std::vector<Cost> result = ComputeShortestPaths(in.source, graph_);
  if (result.size() != static_cast<std::size_t>(in.n)) {
    return false;
  }

  OutType total = 0;
  if (!SumReachableDistances(result, total)) {
    return false;
  }
  GetOutput() = total;
  return true;
}

bool NalitovDDijkstrasAlgorithmTBB::PostProcessingImpl() {
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm
