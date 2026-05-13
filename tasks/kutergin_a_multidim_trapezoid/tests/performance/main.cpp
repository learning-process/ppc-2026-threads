#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"
#include "kutergin_a_multidim_trapezoid/tbb/include/ops_tbb.hpp"
#include "kutergin_a_multidim_trapezoid/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kutergin_a_multidim_trapezoid {

class KuterginATrapezoidPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  static constexpr int kGridSize = 20;
 protected:
  void SetUp() override {
    std::function<double(const std::vector<double>&)> func = [](const std::vector<double> &x) {
      return x[0] + (x[1] * x[1]) + (x[2] * x[2] * x[2]);
    };

    std::vector<std::pair<double, double>> bounds = {{0.0, 2.0}, {0.0, 2.0}, {0.0, 2.0}};

    input_data_ = InType{func, bounds, kGridSize};
    expected_ = 34.6666666667;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output) final {
    constexpr double kTolerance = 1.0;
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

// Настройки для SEQ
const auto kPerfValuesSEQ = ppc::util::TupleToGTestValues(
    ppc::util::MakeAllPerfTasks<InType, KuterginAMultidimTrapezoidSEQ>(
        PPC_SETTINGS_kutergin_a_multidim_trapezoid));

// Настройки для TBB
const auto kPerfValuesTBB = ppc::util::TupleToGTestValues(
    ppc::util::MakeAllPerfTasks<InType, KuterginAMultidimTrapezoidTBB>(
        PPC_SETTINGS_kutergin_a_multidim_trapezoid));

// Регистрация наборов
INSTANTIATE_TEST_SUITE_P(KuterginATrapezoidPerfSEQ, KuterginATrapezoidPerfTest,
                         kPerfValuesSEQ);

INSTANTIATE_TEST_SUITE_P(KuterginATrapezoidPerfTBB, KuterginATrapezoidPerfTest,
                         kPerfValuesTBB);

}  // namespace

}  // namespace kutergin_a_multidim_trapezoid
