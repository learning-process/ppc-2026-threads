#include <gtest/gtest.h>

#include "karpich_i_bitwise_batcher_seq/common/include/common.hpp"
#include "karpich_i_bitwise_batcher_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace karpich_i_bitwise_batcher_seq {

class KarpichIBitwiseBatcherPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1500000;
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

TEST_P(KarpichIBitwiseBatcherPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KarpichIBitwiseBatcherSEQ>(PPC_SETTINGS_karpich_i_bitwise_batcher_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KarpichIBitwiseBatcherPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KarpichIBitwiseBatcherPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace karpich_i_bitwise_batcher_seq
