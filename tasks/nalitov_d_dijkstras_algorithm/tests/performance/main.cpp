#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "nalitov_d_dijkstras_algorithm/all/include/ops_all.hpp"
#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"
#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"
#include "nalitov_d_dijkstras_algorithm/stl/include/ops_stl.hpp"
#include "nalitov_d_dijkstras_algorithm/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

InType MakeRandomSparseGraph(int n, int degree_per_vertex, int seed) {
  InType g;
  g.n = n;
  g.source = 0;
  g.arcs.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(degree_per_vertex));

  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> weight_dist(1, 100);
  std::uniform_int_distribution<int> to_dist(0, n - 1);

  for (int from = 0; from < n; ++from) {
    int added = 0;
    while (added < degree_per_vertex) {
      int to = to_dist(rng);
      if (to != from) {
        g.arcs.push_back(Arc{.from = from, .to = to, .weight = weight_dist(rng)});
        ++added;
      }
    }
  }
  return g;
}

int FindNextVertex(int n, const std::vector<Cost> &dist, const std::vector<char> &visited) {
  int u = -1;
  Cost best = kInf;
  for (int vi = 0; vi < n; ++vi) {
    if (visited[static_cast<std::size_t>(vi)] == 0 && dist[static_cast<std::size_t>(vi)] < best) {
      best = dist[static_cast<std::size_t>(vi)];
      u = vi;
    }
  }
  return u;
}

void RelaxEdges(int u, const std::vector<std::vector<std::pair<int, int>>> &adj, std::vector<Cost> &dist,
                const std::vector<char> &visited) {
  for (const auto &[v, w] : adj[static_cast<std::size_t>(u)]) {
    if (visited[static_cast<std::size_t>(v)] == 0 && dist[static_cast<std::size_t>(u)] <= kInf - w &&
        dist[static_cast<std::size_t>(u)] + w < dist[static_cast<std::size_t>(v)]) {
      dist[static_cast<std::size_t>(v)] = dist[static_cast<std::size_t>(u)] + w;
    }
  }
}

OutType ComputeExpectedOutput(const InType &g) {
  std::vector<Cost> dist(static_cast<std::size_t>(g.n), kInf);
  std::vector<char> visited(static_cast<std::size_t>(g.n), 0);
  dist[static_cast<std::size_t>(g.source)] = 0;

  std::vector<std::vector<std::pair<int, int>>> adj(static_cast<std::size_t>(g.n));
  for (const Arc &a : g.arcs) {
    adj[static_cast<std::size_t>(a.from)].emplace_back(a.to, a.weight);
  }

  for (int step = 0; step < g.n; ++step) {
    int u = FindNextVertex(g.n, dist, visited);
    if (u == -1) {
      break;
    }
    visited[static_cast<std::size_t>(u)] = 1;
    RelaxEdges(u, adj, dist, visited);
  }

  std::int64_t acc = 0;
  for (Cost d : dist) {
    if (d != kInf) {
      if (acc > std::numeric_limits<std::int64_t>::max() - d) {
        return 0;
      }
      acc += d;
    }
  }
  return static_cast<OutType>(acc);
}

constexpr int kGraphSize = 10000;
constexpr int kDegree = 20;
constexpr int kSeed = 42;

}  // namespace

class NalitovDDijkstrasAlgorithmPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  InType input_data{};
  OutType expected_output{};

  void SetUp() override {
    static const InType kStaticInput = MakeRandomSparseGraph(kGraphSize, kDegree, kSeed);
    static const OutType kStaticOutput = ComputeExpectedOutput(kStaticInput);

    input_data = kStaticInput;
    expected_output = kStaticOutput;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output == output_data;
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(NalitovDDijkstrasAlgorithmPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kSeqPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmSeq>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kOmpPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmOmp>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kTbbPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmTBB>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kStlPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmSTL>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmALL>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kPerfTasks = std::tuple_cat(kSeqPerfTasks, kOmpPerfTasks, kTbbPerfTasks, kStlPerfTasks, kAllPerfTasks);

const auto kGtestValues = ppc::util::TupleToGTestValues(kPerfTasks);
const auto kPerfTestName = NalitovDDijkstrasAlgorithmPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NalitovDDijkstrasAlgorithmPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm
