#include <gtest/gtest.h>

#include "dergachev_a_graham_scan/all/include/ops_all.hpp"
#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "dergachev_a_graham_scan/omp/include/ops_omp.hpp"
#include "dergachev_a_graham_scan/seq/include/ops_seq.hpp"
#include "dergachev_a_graham_scan/stl/include/ops_stl.hpp"
#include "dergachev_a_graham_scan/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace dergachev_a_graham_scan {

class DergachevAGrahamScanPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
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

TEST_P(DergachevAGrahamScanPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NesterovATestTaskALL, NesterovATestTaskOMP, DergachevAGrahamScanSEQ,
                                NesterovATestTaskSTL, NesterovATestTaskTBB>(PPC_SETTINGS_dergachev_a_graham_scan);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DergachevAGrahamScanPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DergachevAGrahamScanPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dergachev_a_graham_scan
