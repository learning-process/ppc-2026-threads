#include <gtest/gtest.h>

#include "kamaletdinov_r_bitwise_int_seq/common/include/common.hpp"
#include "kamaletdinov_r_bitwise_int_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kamaletdinov_r_bitwise_int_seq {

class KamaletdinovRBitwiseIntPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1000000;
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

TEST_P(KamaletdinovRBitwiseIntPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KamaletdinovRBitwiseIntSEQ>(PPC_SETTINGS_kamaletdinov_r_bitwise_int_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KamaletdinovRBitwiseIntPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KamaletdinovRBitwiseIntPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamaletdinov_r_bitwise_int_seq
