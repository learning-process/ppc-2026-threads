#include <gtest/gtest.h>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "gasenin_l_djstra/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace gasenin_l_djstra {

class DjkstraPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};
  OutType expected_output_{};

  void SetUp() override {
    input_data_ = 200;
    expected_output_ = input_data_ * (input_data_ - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(DjkstraPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GaseninLDjstraSEQ>(PPC_SETTINGS_gasenin_l_djstra);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = DjkstraPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(DjkstraSeqPerf, DjkstraPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace gasenin_l_djstra
