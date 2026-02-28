#include <gtest/gtest.h>

#include "peterson_r_graham_scan_omp/common/include/common.hpp"
#include "peterson_r_graham_scan_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace peterson_r_graham_scan_omp {

class PetersonRGrahamScanPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 500000;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PetersonRGrahamScanPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PetersonRGrahamScanOMP>(PPC_SETTINGS_peterson_r_graham_scan_omp);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PetersonRGrahamScanPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PetersonRGrahamScanPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace peterson_r_graham_scan_omp
