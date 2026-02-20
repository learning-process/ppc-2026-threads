#include <gtest/gtest.h>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace levonychev_i_radix_batcher_sort {

class LevonychevIRadixBatcherSortRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() > 0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(LevonychevIRadixBatcherSortRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LevonychevIRadixBatcherSortSEQ>(PPC_SETTINGS_levonychev_i_radix_batcher_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LevonychevIRadixBatcherSortRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, LevonychevIRadixBatcherSortRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace levonychev_i_radix_batcher_sort
