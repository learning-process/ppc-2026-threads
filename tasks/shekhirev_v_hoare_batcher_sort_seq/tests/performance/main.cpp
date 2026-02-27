#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "common/include/common.hpp"
#include "seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

class ShekhirevVRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  ShekhirevVRunPerfTest() : kArraySize_(200000) {}

 protected:
  const size_t kArraySize_;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(kArraySize_);
    for (size_t i = 0; i < kArraySize_; ++i) {
      input_data_[i] = static_cast<int>(kArraySize_ - i);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ShekhirevVRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShekhirevHoareBatcherSortSEQ>(PPC_SETTINGS_shekhirev_v_hoare_batcher_sort_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShekhirevVRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShekhirevVRunPerfTest, kGtestValues, kPerfTestName);
}  // namespace

}  // namespace shekhirev_v_hoare_batcher_sort_seq
