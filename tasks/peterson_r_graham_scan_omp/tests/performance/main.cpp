#include <gtest/gtest.h>
#include <omp.h>

#include "peterson_r_graham_scan_omp/common/include/common.hpp"
#include "peterson_r_graham_scan_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace peterson_r_graham_scan_omp {

class PetersonGrahamScannerOMPPerfTests : public ppc::util::BaseRunPerfTests<InputValue, OutputValue> {
  const int kSampleSize_ = 500000;
  InputValue test_input_{};

  void SetUp() override {
    test_input_ = kSampleSize_;
  }

  bool CheckTestOutputData(OutputValue &output_data) final {
    return test_input_ == output_data;
  }

  InputValue GetTestInputData() final {
    return test_input_;
  }
};

TEST_P(PetersonGrahamScannerOMPPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerOMP>(PPC_SETTINGS_peterson_r_graham_scan_omp);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kTestNameGenerator = PetersonGrahamScannerOMPPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, PetersonGrahamScannerOMPPerfTests, kGtestValues, kTestNameGenerator);

}  // namespace

}  // namespace peterson_r_graham_scan_omp
