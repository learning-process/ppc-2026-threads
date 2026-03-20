#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace pankov_a_path_dejikstra {

using Weight = int;
using Vertex = std::size_t;
using Edge = std::tuple<Vertex, Vertex, Weight>;

struct GraphInput {
  Vertex n{};
  Vertex source{};
  std::vector<Edge> edges;
};

constexpr Weight kInfinity = std::numeric_limits<Weight>::max() / 4;

using InType = GraphInput;
using OutType = std::vector<Weight>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace pankov_a_path_dejikstra
