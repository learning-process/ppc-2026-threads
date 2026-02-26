#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "shilin_n_monte_carlo_integration/common/include/common.hpp"
#include "shilin_n_monte_carlo_integration/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shilin_n_monte_carlo_integration {

class ShilinNRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    input_data_ = std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, 10000000,
                                  FuncType::kSumSquares);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &[lower, upper, n, func_type] = input_data_;
    double expected = IntegrandFunction::AnalyticalIntegral(func_type, lower, upper);

    double volume = 1.0;
    for (size_t i = 0; i < lower.size(); ++i) {
      volume *= (upper[i] - lower[i]);
    }
    double epsilon = std::max(volume * 10.0 / std::sqrt(static_cast<double>(n)), 1e-2);
    return std::abs(output_data - expected) <= epsilon;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ShilinNRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShilinNMonteCarloIntegrationSEQ>(PPC_SETTINGS_shilin_n_monte_carlo_integration);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShilinNRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShilinNRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shilin_n_monte_carlo_integration
