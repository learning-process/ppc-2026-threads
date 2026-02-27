#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "shkrebko_m_calc_of_integral_rect/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

class ShkrebkoMRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  static constexpr int kN = 10000;

 protected:
  void SetUp() override {
    input_data_ =
        InType{{{0.0, 1.0}, {0.0, 1.0}}, kN, [](const std::vector<double> &x) { return x[0] * x[0] + x[1] * x[1]; }};
    expected_ = 2.0 / 3.0;
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

TEST_P(ShkrebkoMRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShkrebkoMCalcOfIntegralRectSEQ>(PPC_SETTINGS_shkrebko_m_calc_of_integral_rect);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = ShkrebkoMRunPerfTest::CustomPerfTestName;
INSTANTIATE_TEST_SUITE_P(RunModeTests, ShkrebkoMRunPerfTest, kGtestValues, kPerfTestName);
}  // namespace

}  // namespace shkrebko_m_calc_of_integral_rect
