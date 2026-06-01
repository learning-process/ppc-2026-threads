#include <gtest/gtest.h>

#include <tuple>

#include "peterson_r_graham_scan/all/include/ops_all.hpp"
#include "peterson_r_graham_scan/common/include/common.hpp"
#include "peterson_r_graham_scan/omp/include/ops_omp.hpp"
#include "peterson_r_graham_scan/seq/include/ops_seq.hpp"
#include "peterson_r_graham_scan/stl/include/ops_stl.hpp"
#include "peterson_r_graham_scan/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace peterson_r_graham_scan {

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

const auto kAllPerfTasks = std::tuple_cat(
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerSeq>(PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerOmp>(PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerTbb>(PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerStl>(PPC_SETTINGS_peterson_r_graham_scan),
    ppc::util::MakeAllPerfTasks<InputValue, PetersonGrahamScannerAll>(PPC_SETTINGS_peterson_r_graham_scan));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kTestNameGenerator = PetersonGrahamScannerPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, PetersonGrahamScannerPerfTests, kGtestValues, kTestNameGenerator);

}  // namespace

}  // namespace peterson_r_graham_scan
