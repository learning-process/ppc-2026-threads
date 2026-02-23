#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"
#include "sannikov_i_integrals_rectangle_method/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sannikov_i_integrals_rectangle_method {

class SannikovIRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  static constexpr int kN = 10000;

 protected:
  void SetUp() override {
    input_data_ = InType{[](const std::vector<double> &x) { return std::sin(x[0]) * std::sin(x[1]); },
                         std::vector<std::pair<double, double>>{{0.0, std::numbers::pi}, {0.0, std::numbers::pi}}, kN};
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

TEST_P(SannikovIRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SannikovIIntegralsRectangleMethodSEQ>(
    PPC_SETTINGS_sannikov_i_integrals_rectangle_method);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SannikovIRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SannikovIRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sannikov_i_integrals_rectangle_method
