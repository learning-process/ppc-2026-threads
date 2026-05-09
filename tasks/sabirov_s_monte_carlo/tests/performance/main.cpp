#include <gtest/gtest.h>

#include <cmath>
#include <tuple>

#include "sabirov_s_monte_carlo/common/include/common.hpp"
#include "sabirov_s_monte_carlo/omp/include/ops_seq.hpp"
#include "sabirov_s_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sabirov_s_monte_carlo {

class SabirovSMonteCarloPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    input_data_ = MCInput{
        .lower = {0.0, 0.0, 0.0}, .upper = {1.0, 1.0, 1.0}, .num_samples = 10000000, .func_type = FuncType::kLinear};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    constexpr double kExpected = 1.5;
    double tol = 5.0 / std::sqrt(static_cast<double>(input_data_.num_samples));
    return std::abs(output_data - kExpected) <= tol;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SabirovSMonteCarloPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    std::tuple_cat(ppc::util::MakePerfTaskTuples<SabirovSMonteCarloSEQ, InType>(PPC_SETTINGS_sabirov_s_monte_carlo),
                   ppc::util::MakePerfTaskTuples<SabirovSMonteCarloOMP, InType>(PPC_SETTINGS_sabirov_s_monte_carlo));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SabirovSMonteCarloPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SabirovSMonteCarloPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabirov_s_monte_carlo
