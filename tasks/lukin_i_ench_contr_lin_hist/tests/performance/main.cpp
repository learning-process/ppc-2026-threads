#include <gtest/gtest.h>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "lukin_i_ench_contr_lin_hist/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace lukin_i_ench_contr_lin_hist {

class LukinIPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int image_size_ = 4096;
  InType input_data_{};

  void SetUp() override {
    int count = std::pow(image_size_, 2);

    input_data_.resize(count);
    for (int i = 0; i < count; i++) {
      input_data_[i] = 80 + i % 81;  //[80, 160] - обычное фото
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    unsigned char min = 255;
    unsigned char max = 0;

    for (const auto &elem : input_data_) {
      if (elem < min) {
        min = elem;
      }
      if (elem > max) {
        max = elem;
      }
    }

    float scale = 255.0f / (max - min);
    for (int i = 0; i < static_cast<int>(input_data_.size()); i++) {
      unsigned char expected_value = (input_data_[i] - min) * scale;
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

TEST_P(LukinIPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LukinITestTaskSEQ>(PPC_SETTINGS_lukin_i_ench_contr_lin_hist);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LukinIPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, LukinIPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace lukin_i_ench_contr_lin_hist
