#include <gtest/gtest.h>

#include <cstddef>
#include <utility>

#include "pankov_a_path_dejikstra/common/include/common.hpp"
#include "pankov_a_path_dejikstra/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace pankov_a_path_dejikstra {

class PankovAPathDejikstraRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static InType MakeSparseGraph(std::size_t n) {
    InType graph;
    graph.n = n;
    graph.source = 0;
    graph.edges.reserve((n > 1 ? (n - 1) : 0) + (n > 2 ? (n - 2) : 0));

    for (std::size_t i = 0; i + 1 < n; i++) {
      graph.edges.emplace_back(i, i + 1, 1);
    }
    for (std::size_t i = 0; i + 2 < n; i++) {
      graph.edges.emplace_back(i, i + 2, 2);
    }
    return graph;
  }

  const std::size_t kCount_ = 3000;
  InType input_data_{};

  void SetUp() override {
    input_data_ = MakeSparseGraph(kCount_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.n || output_data[input_data_.source] != 0) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); i++) {
      if (output_data[i] < 0 || std::cmp_not_equal(output_data[i], i)) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PankovAPathDejikstraRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PankovAPathDejikstraSEQ>(PPC_SETTINGS_pankov_a_path_dejikstra);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PankovAPathDejikstraRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PankovAPathDejikstraRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pankov_a_path_dejikstra
