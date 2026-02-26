#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "gonozov_l_bitwise_sorting_double_Batcher_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace gonozov_l_bitwise_sorting_double_Batcher_merge {

class GonozovLBitSortBatcherMergePerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 1000000;
  InType input_data_;

  void SetUp() override {
    std::vector<double> forming_data;
    for (size_t i = 0; i < kCount; i++) {
      forming_data.push_back((i + 120358361) * (i + 1) / (234134) % 1000);
    }
    input_data_ = forming_data;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    for (size_t i = 1; i < kCount; i++) {
      if (output_data[i] < output_data[i - 1]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(GonozovLBitSortBatcherMergePerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GonozovLBitSortBatcherMergeSEQ>(
    PPC_SETTINGS_gonozov_l_bitwise_sorting_double_Batcher_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GonozovLBitSortBatcherMergePerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GonozovLBitSortBatcherMergePerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gonozov_l_bitwise_sorting_double_Batcher_merge
