#include <gtest/gtest.h>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
// #include "perepelkin_i_convex_hull_graham_scan/mpi/include/ops_mpi.hpp"
#include "perepelkin_i_convex_hull_graham_scan/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 100;
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

TEST_P(PerepelkinIConvexHullGrahamScanPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PerepelkinIConvexHullGrahamScanSEQ>(PPC_SETTINGS_perepelkin_i_convex_hull_graham_scan);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PerepelkinIConvexHullGrahamScanPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PerepelkinIConvexHullGrahamScanPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace perepelkin_i_convex_hull_graham_scan
