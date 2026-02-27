#include <gtest/gtest.h>

#include "telnov_a_integral_rectangle/common/include/common.hpp"
#include "telnov_a_integral_rectangle/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace telnov_a_integral_rectangle {

class TelnovAIntegralRectanglePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    input_data_ = {50, 4};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double expected = 4.0 / 2.0;
    return std::abs(output_data - expected) < 1e-2;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(TelnovAIntegralRectanglePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TelnovAIntegralRectangleSEQ>(PPC_SETTINGS_telnov_a_integral_rectangle);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TelnovAIntegralRectanglePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TelnovAIntegralRectanglePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace telnov_a_integral_rectangle
