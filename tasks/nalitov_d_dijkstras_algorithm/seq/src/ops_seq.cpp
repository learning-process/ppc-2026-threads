#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;

struct FrontierNode {
  Cost dist{};
  NodeId node{};
  friend bool operator<(const FrontierNode &a, const FrontierNode &b) {
    if (a.dist != b.dist) {
      return a.dist < b.dist;
    }
    return a.node < b.node;
  }
  friend bool operator>(const FrontierNode &a, const FrontierNode &b) {
    return b < a;
  }
};

bool CheckedSum(std::int64_t acc, Cost addend, std::int64_t &out) {
  const auto x = static_cast<std::int64_t>(addend);
  if (x > 0 && acc > std::numeric_limits<std::int64_t>::max() - x) {
    return false;
  }
  if (x < 0 && acc < std::numeric_limits<std::int64_t>::min() - x) {
    return false;
  }
  out = acc + x;
  return true;
}

std::vector<Cost> FindShortestPaths(NodeId start, const OutgoingTable &graph) {
  const auto vertex_count = graph.size();
  std::vector<Cost> best(vertex_count, kInf);
  std::priority_queue<FrontierNode, std::vector<FrontierNode>, std::greater<>> pq;

  best[static_cast<std::size_t>(start)] = 0;
  pq.push({0, start});

  while (!pq.empty()) {
    const FrontierNode top = pq.top();
    pq.pop();
    const Cost dist_u = top.dist;
    const NodeId u = top.node;
    const auto ui = static_cast<std::size_t>(u);

    if (dist_u != best[ui]) {
      continue;
    }

    for (const auto &neighbor : graph[ui]) {
      const NodeId v = neighbor.first;
      const Cost w = neighbor.second;
      const auto vi = static_cast<std::size_t>(v);

      if (dist_u <= kInf - w && dist_u + w < best[vi]) {
        best[vi] = dist_u + w;
        pq.push({best[vi], v});
      }
    }
  }

  return best;
}

bool AccumulateFiniteDistances(const std::vector<Cost> &best, OutType &sum) {
  std::int64_t acc = 0;
  for (Cost d : best) {
    if (d == kInf) {
      continue;
    }
    if (!CheckedSum(acc, d, acc)) {
      return false;
    }
  }
  if (acc < 0 || acc > std::numeric_limits<OutType>::max()) {
    return false;
  }
  sum = acc;
  return true;
}

}  // namespace

NalitovDDijkstrasAlgorithmSeq::NalitovDDijkstrasAlgorithmSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool NalitovDDijkstrasAlgorithmSeq::ValidationImpl() {
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

  const auto arc_ok = [&in](const Arc &a) {
    return a.from >= 0 && a.to >= 0 && a.from < in.n && a.to < in.n && a.weight >= 0;
  };
  return std::ranges::all_of(in.arcs, arc_ok);
}

bool NalitovDDijkstrasAlgorithmSeq::PreProcessingImpl() {
  const InType &in = GetInput();
  graph_.assign(static_cast<std::size_t>(in.n), {});

  for (const Arc &a : in.arcs) {
    graph_[static_cast<std::size_t>(a.from)].emplace_back(a.to, a.weight);
  }
  GetOutput() = 0;
  return true;
}

bool NalitovDDijkstrasAlgorithmSeq::RunImpl() {
  const InType &in = GetInput();
  if (graph_.size() != static_cast<std::size_t>(in.n)) {
    return false;
  }

  const std::vector<Cost> best = FindShortestPaths(in.source, graph_);
  if (best.size() != static_cast<std::size_t>(in.n)) {
    return false;
  }

  OutType total = 0;
  if (!AccumulateFiniteDistances(best, total)) {
    return false;
  }
  GetOutput() = total;
  return true;
}

bool NalitovDDijkstrasAlgorithmSeq::PostProcessingImpl() {
  return GetOutput() >= 0;
}

}  // namespace nalitov_d_dijkstras_algorithm
