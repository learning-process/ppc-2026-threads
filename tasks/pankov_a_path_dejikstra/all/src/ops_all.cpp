#include "pankov_a_path_dejikstra/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
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

std::vector<std::int64_t> PackEdges(const std::vector<Edge> &edges) {
  std::vector<std::int64_t> packed(edges.size() * 3);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto &[from, to, weight] = edges[i];
    const std::size_t offset = i * 3;
    packed[offset] = static_cast<std::int64_t>(from);
    packed[offset + 1] = static_cast<std::int64_t>(to);
    packed[offset + 2] = static_cast<std::int64_t>(weight);
  }
  return packed;
}

std::vector<Edge> UnpackEdges(const std::vector<std::int64_t> &packed) {
  std::vector<Edge> edges(packed.size() / 3);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const std::size_t offset = i * 3;
    edges[i] = {static_cast<Vertex>(packed[offset]), static_cast<Vertex>(packed[offset + 1]),
                static_cast<Weight>(packed[offset + 2])};
  }
  return edges;
}

AdjList BuildAdjacency(Vertex vertices_count, const std::vector<Edge> &edges) {
  AdjList adjacency(vertices_count);
  for (const auto &[from, to, weight] : edges) {
    adjacency[from].emplace_back(to, weight);
  }
  return adjacency;
}

}  // namespace

PankovAPathDejikstraALL::PankovAPathDejikstraALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool PankovAPathDejikstraALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int is_valid = 1;
  if (rank == 0) {
    const InType &input = GetInput();
    if (input.n == 0 || input.source >= input.n) {
      is_valid = 0;
    } else {
      const auto edge_valid = [&input](const Edge &e) {
        const auto [from, to, weight] = e;
        return from < input.n && to < input.n && weight >= 0;
      };
      is_valid = std::ranges::all_of(input.edges, edge_valid) ? 1 : 0;
    }
  }

  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid == 1;
}

bool PankovAPathDejikstraALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<std::uint64_t, 3> metadata{};
  std::vector<std::int64_t> packed_edges;
  if (rank == 0) {
    const InType &input = GetInput();
    metadata[0] = static_cast<std::uint64_t>(input.n);
    metadata[1] = static_cast<std::uint64_t>(input.source);
    metadata[2] = static_cast<std::uint64_t>(input.edges.size());
    packed_edges = PackEdges(input.edges);
  }

  MPI_Bcast(metadata.data(), static_cast<int>(metadata.size()), MPI_UINT64_T, 0, MPI_COMM_WORLD);

  vertices_count_ = static_cast<Vertex>(metadata[0]);
  source_ = static_cast<Vertex>(metadata[1]);
  packed_edges.resize(static_cast<std::size_t>(metadata[2]) * 3);
  if (!packed_edges.empty()) {
    MPI_Bcast(packed_edges.data(), static_cast<int>(packed_edges.size()), MPI_LONG, 0, MPI_COMM_WORLD);
  }

  adjacency_ = BuildAdjacency(vertices_count_, UnpackEdges(packed_edges));
  GetOutput().clear();
  return true;
}

bool PankovAPathDejikstraALL::RunImpl() {
  OutType distances;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    if (adjacency_.size() != vertices_count_) {
      return false;
    }
    distances = DijkstraSeq(source_, adjacency_);
  }

  auto output_size = static_cast<std::uint64_t>(distances.size());
  MPI_Bcast(&output_size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  GetOutput().assign(static_cast<std::size_t>(output_size), kInfinity);
  if (rank == 0) {
    GetOutput() = std::move(distances);
  }
  if (!GetOutput().empty()) {
    MPI_Bcast(GetOutput().data(), static_cast<int>(GetOutput().size()), MPI_INT, 0, MPI_COMM_WORLD);
  }

  return GetOutput().size() == vertices_count_;
}

bool PankovAPathDejikstraALL::PostProcessingImpl() {
  return GetOutput().size() == vertices_count_;
}

}  // namespace pankov_a_path_dejikstra
