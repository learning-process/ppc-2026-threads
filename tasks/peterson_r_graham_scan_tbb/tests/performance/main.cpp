#include <gtest/gtest.h>

#include "peterson_r_graham_scan_tbb/common/include/common.hpp"
#include "peterson_r_graham_scan_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace peterson_r_graham_scan_tbb {

class PetersonGrahamScannerTBBPerfTests : public ppc::util::BaseRunPerfTests<InputValue, OutputValue> {
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

TEST_P(PetersonGrahamScannerTBBPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerTBB>(PPC_SETTINGS_peterson_r_graham_scan_tbb);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kTestNameGenerator = PetersonGrahamScannerTBBPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, PetersonGrahamScannerTBBPerfTests, kGtestValues, kTestNameGenerator);

}  // namespace

}  // namespace peterson_r_graham_scan_tbb
