#include <gtest/gtest.h>

#include "peterson_r_graham_scan_seq/common/include/common.hpp"
#include "peterson_r_graham_scan_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace peterson_r_graham_scan_seq {

class PetersonGrahamScannerPerfTests : public ppc::util::BaseRunPerfTests<InputValue, OutputValue> {
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

TEST_P(PetersonGrahamScannerPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScanner>(PPC_SETTINGS_peterson_r_graham_scan_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kTestNameGenerator = PetersonGrahamScannerPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, PetersonGrahamScannerPerfTests, kGtestValues, kTestNameGenerator);

}  // namespace

}  // namespace peterson_r_graham_scan_seq
