#include "pankov_a_path_dejikstra/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

#include "pankov_a_path_dejikstra/common/include/common.hpp"
#include "util/include/util.hpp"

namespace pankov_a_path_dejikstra {
namespace {

using AdjList = std::vector<std::vector<std::pair<Vertex, Weight>>>;

OutType DijkstraSeq(Vertex source, const AdjList &adjacency) {
  OutType distance(adjacency.size(), kInfinity);
  using QueueNode = std::pair<Weight, Vertex>;
  std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>> min_queue;

  distance[source] = 0;
  min_queue.emplace(0, source);

  while (!min_queue.empty()) {
    const auto [current_dist, u] = min_queue.top();
    min_queue.pop();

    if (current_dist != distance[u]) {
      continue;
    }

    for (const auto &[v, weight] : adjacency[u]) {
      if (current_dist <= kInfinity - weight && current_dist + weight < distance[v]) {
        distance[v] = current_dist + weight;
        min_queue.emplace(distance[v], v);
      }
    }
  }

  return distance;
}

}  // namespace

PankovAPathDejikstraOMP::PankovAPathDejikstraOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PankovAPathDejikstraOMP::ValidationImpl() {
  const InType &input = GetInput();
  if (input.n == 0 || input.source >= input.n) {
    return false;
  }

  const auto edge_valid = [&input](const Edge &e) {
    const auto [from, to, weight] = e;
    return from < input.n && to < input.n && weight >= 0;
  };
  return std::ranges::all_of(input.edges, edge_valid);
}

bool PankovAPathDejikstraOMP::PreProcessingImpl() {
  const InType &input = GetInput();
  adjacency_.assign(input.n, {});
  auto &adjacency = adjacency_;
  const auto vertices_count = static_cast<std::ptrdiff_t>(input.n);
  const int omp_threads = ppc::util::GetNumThreads();
  if (omp_threads <= 0) {
    return false;
  }

  std::vector<Edge> sorted_edges = input.edges;
  std::ranges::sort(sorted_edges, [](const Edge &lhs, const Edge &rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });

  std::vector<std::size_t> offsets(input.n + 1, 0);
  std::size_t edge_idx = 0;
  for (Vertex from = 0; from < input.n; from++) {
    offsets[from] = edge_idx;
    while (edge_idx < sorted_edges.size() && std::get<0>(sorted_edges[edge_idx]) == from) {
      edge_idx++;
    }
  }
  offsets[input.n] = edge_idx;

#pragma omp parallel for default(none) schedule(static) num_threads(omp_threads) \
    shared(adjacency, offsets, sorted_edges, vertices_count)
  for (std::ptrdiff_t from = 0; from < vertices_count; ++from) {
    const std::size_t begin = offsets[static_cast<std::size_t>(from)];
    const std::size_t end = offsets[static_cast<std::size_t>(from) + 1];
    adjacency[static_cast<std::size_t>(from)].reserve(end - begin);
    for (std::size_t i = begin; i < end; i++) {
      const auto &[edge_from, to, weight] = sorted_edges[i];
      (void)edge_from;
      adjacency[static_cast<std::size_t>(from)].emplace_back(to, weight);
    }
  }

  GetOutput().clear();
  return true;
}

bool PankovAPathDejikstraOMP::RunImpl() {
  const InType &input = GetInput();
  if (adjacency_.size() != input.n) {
    return false;
  }
  GetOutput() = DijkstraSeq(input.source, adjacency_);
  return GetOutput().size() == input.n;
}

bool PankovAPathDejikstraOMP::PostProcessingImpl() {
  return GetOutput().size() == GetInput().n;
}

}  // namespace pankov_a_path_dejikstra
