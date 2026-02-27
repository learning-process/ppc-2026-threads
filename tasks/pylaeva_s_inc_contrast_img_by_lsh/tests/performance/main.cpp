#include <gtest/gtest.h>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_;

  void SetUp() override {

    input_data_.resize(kCount_);
    for (int i = 0; i < kCount_; i++) {
      input_data_[i] = 80 + (i % 81);  //[80, 160] - обычное фото
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto min_it = std::ranges::min_element(input_data_.begin(), input_data_.end());
    auto max_it = std::ranges::max_element(input_data_.begin(), input_data_.end());

    unsigned char min = *min_it;
    unsigned char max = *max_it;

    float scale = 255.0F / static_cast<float>(max - min);

    int size = static_cast<int>(input_data_.size());

    for (int i = 0; i < size; i++) {
      auto expected_value = static_cast<unsigned char>(static_cast<float>(input_data_[i] - min) * scale);
      if (output_data[i] != expected_value) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PylaevaSRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PylaevaSIncContrastImgByLshSEQ>(PPC_SETTINGS_pylaeva_s_inc_contrast_img_by_lsh);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PylaevaSRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PylaevaSRunPerfTests, PylaevaSRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
