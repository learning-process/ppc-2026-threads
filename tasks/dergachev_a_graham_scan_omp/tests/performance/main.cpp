#include <gtest/gtest.h>

#include "dergachev_a_graham_scan_omp/common/include/common.hpp"
#include "dergachev_a_graham_scan_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace dergachev_a_graham_scan_omp {

class DergachevAGrahamScanPerfTestsOmp : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(DergachevAGrahamScanPerfTestsOmp, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DergachevAGrahamScanOMP>(PPC_SETTINGS_dergachev_a_graham_scan_omp);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DergachevAGrahamScanPerfTestsOmp::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DergachevAGrahamScanPerfTestsOmp, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dergachev_a_graham_scan_omp
