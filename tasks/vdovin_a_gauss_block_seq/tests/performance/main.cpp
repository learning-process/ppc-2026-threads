#include <gtest/gtest.h>

#include "util/include/perf_test_util.hpp"
#include "vdovin_a_gauss_block_seq/common/include/common.hpp"
#include "vdovin_a_gauss_block_seq/seq/include/ops_seq.hpp"

namespace vdovin_a_gauss_block_seq {

class VdovinAGaussBlockPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1500;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == 100;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(VdovinAGaussBlockPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, VdovinAGaussBlockSEQ>(PPC_SETTINGS_vdovin_a_gauss_block_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = VdovinAGaussBlockPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, VdovinAGaussBlockPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace vdovin_a_gauss_block_seq
