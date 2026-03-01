#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"
#include "eremin_v_integrals_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace eremin_v_integrals_monte_carlo {

class EreminVRunPerfTestsThreadsIntegralsMonteCarlo : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    MonteCarloInput input;
    input.bounds = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
    input.samples = 5'000'000;
    input.func = [](const std::vector<double> &x) { return x[0] + x[1] + x[2]; };

    input_data_ = input;

    expected_result_ = 1.5;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    double tolerance = 1e-2;
    return std::abs(output_data - expected_result_) <= tolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_{};
};

TEST_P(EreminVRunPerfTestsThreadsIntegralsMonteCarlo, IntegralsMonteCarloPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, EreminVIntegralsMonteCarloSEQ>(PPC_SETTINGS_eremin_v_integrals_monte_carlo);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = EreminVRunPerfTestsThreadsIntegralsMonteCarlo::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(IntegralsMonteCarloTestsPerf, EreminVRunPerfTestsThreadsIntegralsMonteCarlo, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace eremin_v_integrals_monte_carlo
