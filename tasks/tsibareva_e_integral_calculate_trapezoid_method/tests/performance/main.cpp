#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaERunPerfTestThreads : public ppc::util::BaseRunPerfTests<Integral, double> {
  Integral input_data_{};

  void SetUp() override {
    input_data_.dim = 3;
    input_data_.lo = {0.0, 0.0, 0.0};
    input_data_.hi = {1.0, 1.0, 1.0};
    input_data_.steps = {300, 300, 300};
    input_data_.f = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); };
  }

  bool CheckTestOutputData(double &output_data) final {
    double expected = 1.0;
    return std::fabs(output_data - expected) < 1e-4;
  }

  Integral GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TsibarevaERunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<Integral, TsibarevaEIntegralCalculateTrapezoidMethodSEQ>(
    PPC_SETTINGS_tsibareva_e_integral_calculate_trapezoid_method);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TsibarevaERunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TsibarevaERunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
