#include <gtest/gtest.h>

#include <cmath>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "lukin_i_ench_contr_lin_hist/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace lukin_i_ench_contr_lin_hist {

class LukinIPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int image_size_ = 4096;
  InType input_data_;

  void SetUp() override {
    int count = static_cast<int>(std::pow(image_size_, 2));

    input_data_.resize(count);
    for (int i = 0; i < count; i++) {
      input_data_[i] = 80 + (i % 81);  //[80, 160] - обычное фото
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto min_it = std::min(input_data_.begin(), input_data_.end());
    auto max_it = std::max(input_data_.begin(), input_data_.end());

    unsigned char min = *min_it;
    unsigned char max = *max_it;

    float scale = 255.0F / static_cast<float>(max - min);

    int size = static_cast<int>(input_data_.size());

    for (int i = 0; i < size; i++) {
      unsigned char expected_value = static_cast<unsigned char>(static_cast<float>(input_data_[i] - min) * scale);
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
