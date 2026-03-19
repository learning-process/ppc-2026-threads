#include <gtest/gtest.h>

#include <cstdint>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace levonychev_i_radix_batcher_sort {

class LevonychevIRadixBatcherSortRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int a = 1664525;
    const int c = 1013904223;
    const int m = INT32_MAX;
    int seed = 15;
    const int size = 10'000'000;
    input_data_.clear();
    input_data_.resize(size);

    for (int i = 0; i < size; ++i) {
      seed = (a * seed + c) % m;
      input_data_[i] = (seed % 2000) - 1000;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
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
