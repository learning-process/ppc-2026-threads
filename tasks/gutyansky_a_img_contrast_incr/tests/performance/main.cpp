#include <gtest/gtest.h>

#include <cstddef>
#include <limits>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"
#include "gutyansky_a_img_contrast_incr/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace gutyansky_a_img_contrast_incr {

class GutyanskyARunPerfTestsImgContrastIncr : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const size_t kCount_ = 100000000;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(kCount_);

    input_data_[0] = std::numeric_limits<InType::value_type>::max();

    for (size_t i = 1; i < kCount_; i += 2) {
      input_data_[i] = (i - 1) % std::numeric_limits<InType::value_type>::max();
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(GutyanskyARunPerfTestsImgContrastIncr, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, GutyanskyAImgContrastIncrSEQ>(PPC_SETTINGS_gutyansky_a_img_contrast_incr);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GutyanskyARunPerfTestsImgContrastIncr::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GutyanskyARunPerfTestsImgContrastIncr, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gutyansky_a_img_contrast_incr
