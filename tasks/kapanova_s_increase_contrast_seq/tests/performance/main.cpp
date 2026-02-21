#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

#include "kapanova_s_increase_contrast_seq/common/include/common.hpp"
#include "kapanova_s_increase_contrast_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kapanova_s_increase_contrast_seq {

class KapanovaSIncreaseContrastPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 1024;
  InType input_data_{};

  void SetUp() override {
    input_data_.resize(kSize_ * kSize_);
    for (size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = static_cast<uint8_t>(i % 256);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (input_data_.empty() || output_data.empty()) {
      return false;
    }
    if (input_data_.size() != output_data.size()) {
      return false;
    }

    auto min_out = *std::min_element(output_data.begin(), output_data.end());
    auto max_out = *std::max_element(output_data.begin(), output_data.end());

    return (min_out == 0 && max_out == 255);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KapanovaSIncreaseContrastPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KapanovaSIncreaseContrastSEQ>(PPC_SETTINGS_kapanova_s_increase_contrast_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KapanovaSIncreaseContrastPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTests, KapanovaSIncreaseContrastPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kapanova_s_increase_contrast_seq
