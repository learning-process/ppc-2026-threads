#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "bortsova_a_integrals_rectangle_seq/common/include/common.hpp"
#include "bortsova_a_integrals_rectangle_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace bortsova_a_integrals_rectangle_seq {

class BortsovaAIntegralsRectanglePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};
  double expected_ = 0.0;

  void SetUp() override {
    input_data_.func = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); };
    input_data_.lower_bounds = {0.0, 0.0, 0.0};
    input_data_.upper_bounds = {1.0, 1.0, 1.0};
    input_data_.num_steps = 200;
    expected_ = 1.0;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_) < 1e-2;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BortsovaAIntegralsRectanglePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, BortsovaAIntegralsRectangleSEQ>(
    PPC_SETTINGS_bortsova_a_integrals_rectangle_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BortsovaAIntegralsRectanglePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BortsovaAIntegralsRectanglePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace bortsova_a_integrals_rectangle_seq
