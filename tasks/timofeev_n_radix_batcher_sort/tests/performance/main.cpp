#include <gtest/gtest.h>

#include <utility>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"
#include "timofeev_n_radix_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    int size = 10000;
    input_data_.resize(size, 0);
    for (int i = 0; std::cmp_less(i, static_cast<int>(input_data_.size())); i++) {
      input_data_[i] = i * (i % 2 == 0 ? 1 : -1);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    for (int i = 0; std::cmp_less(i, static_cast<int>(output_data.size() - 1)); i++) {
      if (output_data[i] > output_data[i + 1]) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TimofeevRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TimofeevNRadixBatcherSEQ>(PPC_SETTINGS_timofeev_n_radix_batcher_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TimofeevRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TimofeevRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace timofeev_n_radix_batcher_sort_threads
