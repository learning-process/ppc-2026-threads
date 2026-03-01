#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <tuple>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"
#include "kutergin_a_multidim_trapezoid/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kutergin_a_multidim_trapezoid {

class KuterginATrapezoidPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  static constexpr int kGridSize = 800;

 protected:
  void SetUp() override {
    auto heavy_func = [](const std::vector<double> &x) { return ((std::sin(x[0]) * std::sin(x[1])) + std::exp(x[2])); };

    std::vector<std::pair<double, double>> bounds = {{0.0, std::numbers::pi}, {0.0, std::numbers::pi}, {0.0, 1.0}};

    input_data_ = InType{heavy_func, bounds, kGridSize};

    expected_ = ((4.0 * 1.0) + ((std::numbers::pi * std::numbers::pi) * (std::numbers::e - 1.0)));
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output) final {
    constexpr double kTolerance = 1e-2;
    return (std::fabs(output - expected_) <= kTolerance);
  }

 private:
  InType input_data_;
  OutType expected_ = 0.0;
};

TEST_P(KuterginATrapezoidPerfTest, PerformanceModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KuterginAMultidimTrapezoidSEQ>(PPC_SETTINGS_kutergin_a_multidim_trapezoid);

const auto kPerfValues = ppc::util::TupleToGTestValues(kPerfTasks);

INSTANTIATE_TEST_SUITE_P(KuterginATrapezoidPerfSuite, KuterginATrapezoidPerfTest, kPerfValues,
                         KuterginATrapezoidPerfTest::CustomPerfTestName);

}  // namespace

}  // namespace kutergin_a_multidim_trapezoid
