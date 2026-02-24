#include <gtest/gtest.h>

#include <cmath>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"
#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

class KiselevPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  InType input_data;

  void SetUp() override {
    input_data.left_bounds = {0.0, 0.0};
    input_data.right_bounds = {1.0, 1.0};
    input_data.step_n_size = {1000, 1000};
    input_data.type_function = 0;
    input_data.epsilon = 0.0;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::isfinite(output_data);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(KiselevPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KiselevITestTaskSEQ>(
    PPC_SETTINGS_kiselev_i_trapezoidal_method_for_multidimensional_integrals);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KiselevPerfTests::CustomPerfTestName;

namespace {
INSTANTIATE_TEST_SUITE_P(KiselevPerfTests, KiselevPerfTests, kGtestValues, kPerfTestName);
}  // namespace
}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
