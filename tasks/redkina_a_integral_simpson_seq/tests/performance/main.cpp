#include <gtest/gtest.h>

#include <cmath>
#include <utility>
#include <vector>

#include "redkina_a_integral_simpson_seq/common/include/common.hpp"
#include "redkina_a_integral_simpson_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace redkina_a_integral_simpson_seq {

class RedkinaAIntegralSimpsonPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    auto func = [](const std::vector<double> &x) { return std::exp(-((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]))); };
    std::vector<double> a = {-1.0, -1.0, -1.0};
    std::vector<double> b = {1.0, 1.0, 1.0};
    std::vector<int> n = {200, 200, 200};

    input_data_ = InputData{.func = std::move(func), .a = std::move(a), .b = std::move(b), .n = std::move(n)};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::isfinite(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(RedkinaAIntegralSimpsonPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RedkinaAIntegralSimpsonSEQ>(PPC_SETTINGS_redkina_a_integral_simpson_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RedkinaAIntegralSimpsonPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RedkinaAIntegralSimpsonPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace redkina_a_integral_simpson_seq
