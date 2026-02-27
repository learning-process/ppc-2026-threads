#include <gtest/gtest.h>

#include "nalitov_d_dijkstras_algorithm_seq/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nalitov_d_dijkstras_algorithm_seq {

class NalitovDDijkstrasAlgorithmSeqPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(NalitovDDijkstrasAlgorithmSeqPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmSeq>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NalitovDDijkstrasAlgorithmSeqPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NalitovDDijkstrasAlgorithmSeqPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm_seq
