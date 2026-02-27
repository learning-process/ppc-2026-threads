#include <gtest/gtest.h>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "chernov_t_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(ChernovTRadixSortPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ChernovTRadixSortSEQ>(PPC_SETTINGS_chernov_t_radix_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChernovTRadixSortPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ChernovTRadixSortPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chernov_t_radix_sort
