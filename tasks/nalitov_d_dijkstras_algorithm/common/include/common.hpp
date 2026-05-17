#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace nalitov_d_dijkstras_algorithm {

/// Идентификатор вершины (нумерация 0 .. vertex_count-1).
using NodeId = int;
/// Вес ребра и значение кратчайшего пути в очереди — один тип, чтобы не смешивать переполнения.
using Cost = int;

/// Ориентированное ребро с неотрицательным весом.
struct Arc {
  NodeId from{};
  NodeId to{};
  Cost weight{};
};

/// Описание задачи: размер, старт и явный список рёбер.
struct GraphDescriptor {
  int n{};
  int source{};
  std::vector<Arc> arcs;
};

/// «Бесконечность» для недостижимых вершин (запас под сумму в релаксациях).
constexpr Cost kInf = std::numeric_limits<Cost>::max() / 2;

using InType = GraphDescriptor;
using OutType = std::int64_t;
using TestType = std::tuple<InType, OutType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nalitov_d_dijkstras_algorithm
