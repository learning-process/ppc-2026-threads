#include <gtest/gtest.h>

#include <cmath>
#include <utility>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"
#include "redkina_a_integral_simpson/omp/include/ops_omp.hpp"
#include "redkina_a_integral_simpson/seq/include/ops_seq.hpp"
#include "redkina_a_integral_simpson/stl/include/ops_stl.hpp"
#include "redkina_a_integral_simpson/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace redkina_a_integral_simpson {

class RedkinaAIntegralSimpsonPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    auto func = [](const std::vector<double> &x) { return std::exp(-((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]))); };
    std::vector<double> a = {-1.0, -1.0, -1.0};
    std::vector<double> b = {1.0, 1.0, 1.0};
    std::vector<int> n = {100, 100, 100};

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

const auto kAllPerfTasksSeq =
    ppc::util::MakeAllPerfTasks<InType, RedkinaAIntegralSimpsonSEQ>(PPC_SETTINGS_redkina_a_integral_simpson);

const auto kAllPerfTasksOmp =
    ppc::util::MakeAllPerfTasks<InType, RedkinaAIntegralSimpsonOMP>(PPC_SETTINGS_redkina_a_integral_simpson);

const auto kAllPerfTasksTbb =
    ppc::util::MakeAllPerfTasks<InType, RedkinaAIntegralSimpsonTBB>(PPC_SETTINGS_redkina_a_integral_simpson);

const auto kAllPerfTasksStl =
    ppc::util::MakeAllPerfTasks<InType, RedkinaAIntegralSimpsonSTL>(PPC_SETTINGS_redkina_a_integral_simpson);

const auto kGtestValuesSeq = ppc::util::TupleToGTestValues(kAllPerfTasksSeq);
const auto kGtestValuesOmp = ppc::util::TupleToGTestValues(kAllPerfTasksOmp);
const auto kGtestValuesTbb = ppc::util::TupleToGTestValues(kAllPerfTasksTbb);
const auto kGtestValuesStl = ppc::util::TupleToGTestValues(kAllPerfTasksStl);

const auto kPerfTestName = RedkinaAIntegralSimpsonPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTestsSeq, RedkinaAIntegralSimpsonPerfTests, kGtestValuesSeq, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(RunModeTestsOmp, RedkinaAIntegralSimpsonPerfTests, kGtestValuesOmp, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(RunModeTestsTbb, RedkinaAIntegralSimpsonPerfTests, kGtestValuesTbb, kPerfTestName);
INSTANTIATE_TEST_SUITE_P(RunModeTestsStl, RedkinaAIntegralSimpsonPerfTests, kGtestValuesStl, kPerfTestName);

}  // namespace

}  // namespace redkina_a_integral_simpson
