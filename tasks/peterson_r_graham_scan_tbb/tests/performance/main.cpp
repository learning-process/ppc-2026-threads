#include <gtest/gtest.h>

#include <tuple>

#include "peterson_r_graham_scan/common/include/common.hpp"
#include "peterson_r_graham_scan/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace peterson_r_graham_scan {

class PetersonGrahamScannerTbbPerfTests : public ppc::util::BaseRunPerfTests<InputValue, OutputValue> {
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

TEST_P(PetersonGrahamScannerTbbPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerTbb>(PPC_SETTINGS_peterson_r_graham_scan);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kTestNameGenerator = PetersonGrahamScannerTbbPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, PetersonGrahamScannerTbbPerfTests, kGtestValues, kTestNameGenerator);

}  // namespace

}  // namespace peterson_r_graham_scan
