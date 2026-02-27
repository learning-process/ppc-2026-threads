#include <gtest/gtest.h>

#include "dergachev_a_graham_scan_stl/common/include/common.hpp"
#include "dergachev_a_graham_scan_stl/stl/include/ops_stl.hpp"
#include "util/include/perf_test_util.hpp"

namespace dergachev_a_graham_scan_stl {

class DergachevAGrahamScanPerfTestsStl : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(DergachevAGrahamScanPerfTestsStl, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DergachevAGrahamScanSTL>(PPC_SETTINGS_dergachev_a_graham_scan_stl);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DergachevAGrahamScanPerfTestsStl::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DergachevAGrahamScanPerfTestsStl, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dergachev_a_graham_scan_stl
