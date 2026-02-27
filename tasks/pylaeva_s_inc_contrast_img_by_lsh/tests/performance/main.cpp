#include <gtest/gtest.h>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "pylaeva_s_inc_contrast_img_by_lsh/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const size_t kCount_ = 10000000;
  InType input_data_;

  void SetUp() override {

    input_data_.resize(kCount_);
    
    // Создаем градиент: значения от 0 до 255
    for (size_t i = 0; i < kCount_; ++i) {
      input_data_[i] = static_cast<uint8_t>((i * 255) / kCount_);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Для градиента проверяем, что min стал 0, max стал 255
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    
    auto [out_min, out_max] = std::ranges::minmax(output_data);
    return (out_min == 0 && out_max == 255);
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
