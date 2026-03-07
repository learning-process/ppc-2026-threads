#include <gtest/gtest.h>

#include <cstddef>

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"
#include "batushin_i_incr_contrast_with_lhs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace batushin_i_incr_contrast_with_lhs {

class BatushinIRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const size_t kPixelsCount_ = static_cast<size_t>(8192) * 8192;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(kPixelsCount_);
    for (size_t i = 0; i < kPixelsCount_; i++) {
      input_data_[i] = static_cast<unsigned char>(100 + (i % 51));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty() && output_data.size() == kPixelsCount_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BatushinIRunPerfTestThreads, IncreaseContrastTest) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BatushinITestTaskSEQ>(PPC_SETTINGS_batushin_i_incr_contrast_with_lhs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BatushinIRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTest, BatushinIRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace batushin_i_incr_contrast_with_lhs
