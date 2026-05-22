#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

using OutgoingTable = std::vector<std::vector<std::pair<NodeId, Cost>>>;

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

std::int64_t FindNextNode(const std::vector<Cost> &best, const std::vector<char> &visited) {
  const auto vertex_count = best.size();
  int best_c = kInf;
  std::int64_t best_v = -1;

  for (std::size_t i = 0; i < vertex_count; ++i) {
    if (visited[i] == 0 && best[i] < best_c) {
      best_c = best[i];
      best_v = static_cast<std::int64_t>(i);
    }
  }
  return best_v;
}

void RelaxNeighborsSeq(std::size_t u, Cost u_dist, const OutgoingTable &graph, std::vector<Cost> &best,
                       const std::vector<char> &visited) {
  for (const auto &neighbor : graph[u]) {
    const NodeId v = neighbor.first;
    const Cost w = neighbor.second;
    const auto vi = static_cast<std::size_t>(v);

    if (visited[vi] == 0 && u_dist <= kInf - w && u_dist + w < best[vi]) {
      best[vi] = u_dist + w;
    }
  }
}

std::vector<Cost> FindShortestPaths(NodeId start, const OutgoingTable &graph) {
  const auto vertex_count = graph.size();
  std::vector<Cost> best(vertex_count, kInf);
  std::vector<char> visited(vertex_count, 0);

  best[static_cast<std::size_t>(start)] = 0;

  for (std::size_t step = 0; step < vertex_count; ++step) {
    std::int64_t u = FindNextNode(best, visited);
    if (u == -1) {
      break;
    }

    const auto ui = static_cast<std::size_t>(u);
    visited[ui] = 1;

    RelaxNeighborsSeq(ui, best[ui], graph, best, visited);
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
  if (in.n <= 0 || in.n > kMaxVertices || in.source < 0 || in.source >= in.n) {
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
