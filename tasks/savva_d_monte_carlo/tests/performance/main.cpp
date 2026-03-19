#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"
#include "savva_d_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace savva_d_monte_carlo {

class SavvaDRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  // const int kCount_ = 200;
  InType input_data_;

  void SetUp() override {
    std::vector<double> lower_bounds = {-10.0, -10.0, -10.0};
    std::vector<double> upper_bounds = {10.0, 10.0, 10.0};
    uint64_t num_points = 50000000;
    auto f = [](const std::vector<double> &) { return 1.0; };
    input_data_ = InputData(lower_bounds, upper_bounds, num_points, std::move(f));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    double expected = 8000.0;
    double tolerance = 0.1;
    return std::abs(output_data - expected) / expected <= tolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SavvaDRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SavvaDMonteCarloSEQ>(PPC_SETTINGS_savva_d_monte_carlo);
// SavvaDMonteCarloALL, SavvaDMonteCarloOMP,  SavvaDMonteCarloSTL, SavvaDMonteCarloTBB

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SavvaDRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SavvaDRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace savva_d_monte_carlo
