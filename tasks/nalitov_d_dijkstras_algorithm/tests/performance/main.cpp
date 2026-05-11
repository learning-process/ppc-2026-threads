#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <tuple>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"
#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"
#include "nalitov_d_dijkstras_algorithm/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace nalitov_d_dijkstras_algorithm {

namespace {

InType MakeStarFromZero(int n) {
  InType g;
  g.n = n;
  g.source = 0;
  g.arcs.reserve(static_cast<std::size_t>(std::max(0, n - 1)));
  for (int i = 1; i < n; ++i) {
    g.arcs.push_back(Arc{.from = 0, .to = i, .weight = i});
  }
  return g;
}

}  // namespace

class NalitovDDijkstrasAlgorithmPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kGraphSize = 150;
  InType input_data_{};
  OutType expected_output_{};

  void SetUp() override {
    input_data_ = MakeStarFromZero(kGraphSize);
    expected_output_ = static_cast<OutType>(kGraphSize) * (kGraphSize - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
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
const auto kAllPerfTasks = std::tuple_cat(kSeqPerfTasks, kOmpPerfTasks, kTbbPerfTasks);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NalitovDDijkstrasAlgorithmPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NalitovDDijkstrasAlgorithmPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm
