#include <gtest/gtest.h>

#include "gasenin_l_djstra_tbb/common/include/common.hpp"
#include "gasenin_l_djstra_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraTbbPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(GaseninLDjstraTbbPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GaseninLDjstraTBB>(PPC_SETTINGS_gasenin_l_djstra_tbb);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = GaseninLDjstraTbbPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(DjkstraTbbPerf, GaseninLDjstraTbbPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace gasenin_l_djstra
