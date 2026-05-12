#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"
#include "dergynov_s_integrals_multistep_rectangle/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace dergynov_s_integrals_multistep_rectangle {

class DergynovSIntegralsRectanglePerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  static constexpr int kN = 1000;

 protected:
  void SetUp() override {
    auto func = [](const std::vector<double> &x) { return std::sin(x[0]) * std::sin(x[1]); };
    std::vector<std::pair<double, double>> borders = {{0.0, std::numbers::pi}, {0.0, std::numbers::pi}};
    input_data_ = InType{func, borders, kN};
    expected_ = 4.0;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-3;
    return std::fabs(output_data - expected_) <= eps;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_ = 0.0;
};

TEST_P(DergynovSIntegralsRectanglePerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, DergynovSIntegralsMultistepRectangleSEQ>(
    PPC_SETTINGS_dergynov_s_integrals_multistep_rectangle);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DergynovSIntegralsRectanglePerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DergynovSIntegralsRectanglePerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dergynov_s_integrals_multistep_rectangle
