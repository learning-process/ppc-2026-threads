#include <gtest/gtest.h>

#include <cmath>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"
#include "afanasyev_a_integ_rect_method/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace afanasyev_a_integ_rect_method {

class AfanasyevAIntegRectMethodPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  const InType kDefaultN = 200;

 protected:
  InType n_ = kDefaultN;

  void SetUp() override { n_ = kDefaultN; }

  bool CheckTestOutputData(OutType &output_data) final {
    constexpr int kDim = 3;
    const double i1 = 0.5 * std::sqrt(M_PI) * std::erf(1.0);
    const double expected = std::pow(i1, kDim);
    const double tol = 5e-3;
    return std::fabs(output_data - expected) <= tol;
  }

  InType GetTestInputData() final { return n_; }
};

TEST_P(AfanasyevAIntegRectMethodPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, AfanasyevAIntegRectMethodSEQ>(
    PPC_SETTINGS_afanasyev_a_integ_rect_method);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = AfanasyevAIntegRectMethodPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerf, AfanasyevAIntegRectMethodPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace afanasyev_a_integ_rect_method
