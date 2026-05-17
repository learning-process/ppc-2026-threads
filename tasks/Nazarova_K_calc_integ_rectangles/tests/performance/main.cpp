#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Nazarova_K_calc_integ_rectangles/common/include/common.hpp"
#include "Nazarova_K_calc_integ_rectangles/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nazarova_k_calc_integ_rectangles {

class NazarovaKCalcIntegRectanglesRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    input_data_ = InType{Integrand, {0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, {120U, 120U, 120U}};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    constexpr double kExpected = 42.0;
    return std::abs(output_data - kExpected) < 1e-10;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  static double Integrand(const std::vector<double> &point) {
    return point[0] + (2.0 * point[1]) + (3.0 * point[2]);
  }
};

TEST_P(NazarovaKCalcIntegRectanglesRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NazarovaKCalcIntegRectanglesSEQ>(PPC_SETTINGS_Nazarova_K_calc_integ_rectangles);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NazarovaKCalcIntegRectanglesRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NazarovaKCalcIntegRectanglesRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nazarova_k_calc_integ_rectangles
