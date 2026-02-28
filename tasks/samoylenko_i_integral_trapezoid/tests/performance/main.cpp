#include <gtest/gtest.h>

#include <cmath>

#include "samoylenko_i_integral_trapezoid/common/include/common.hpp"
#include "samoylenko_i_integral_trapezoid/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace samoylenko_i_integral_trapezoid {

class SamoylenkoIIntegralTrapezoidPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    input_data_ = InType{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {260, 260, 260}, 2};
  }

  bool CheckTestOutputData(OutType &output_data) override {
    double expected = 1.0;
    return std::abs(output_data - expected) < 1e-4;
  }

  InType GetTestInputData() override {
    return input_data_;
  }
};

TEST_P(SamoylenkoIIntegralTrapezoidPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SamoylenkoIIntegralTrapezoidSEQ>(PPC_SETTINGS_samoylenko_i_integral_trapezoid);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SamoylenkoIIntegralTrapezoidPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SamoylenkoIIntegralTrapezoidPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace samoylenko_i_integral_trapezoid
