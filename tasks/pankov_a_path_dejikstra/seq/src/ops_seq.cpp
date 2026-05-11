#include "pankov_a_path_dejikstra/seq/include/ops_seq.hpp"

#include <algorithm>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

#include "pankov_a_path_dejikstra/common/include/common.hpp"

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

PankovAPathDejikstraSEQ::PankovAPathDejikstraSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PankovAPathDejikstraSEQ::ValidationImpl() {
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

bool PankovAPathDejikstraSEQ::PreProcessingImpl() {
  const InType &input = GetInput();
  adjacency_.assign(input.n, {});
  for (const auto &[from, to, weight] : input.edges) {
    adjacency_[from].emplace_back(to, weight);
  }
  GetOutput().clear();
  return true;
}

bool PankovAPathDejikstraSEQ::RunImpl() {
  const InType &input = GetInput();
  if (adjacency_.size() != input.n) {
    return false;
  }
  GetOutput() = DijkstraSeq(input.source, adjacency_);
  return GetOutput().size() == input.n;
}

bool PankovAPathDejikstraSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().n;
}

}  // namespace pankov_a_path_dejikstra
