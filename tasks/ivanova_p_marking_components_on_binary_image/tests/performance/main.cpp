#include <gtest/gtest.h>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 512;
  InType input_data_;

  void SetUp() override {
    // Создаем тестовое изображение
    test_image.width = kSize_;
    test_image.height = kSize_;
    test_image.data.resize(kSize_ * kSize_);

    for (int y = 0; y < kSize_; ++y) {
      for (int x = 0; x < kSize_; ++x) {
        int idx = y * kSize_ + x;
        uint8_t pixel = 0;

        if ((x > 50 && x < 150 && y > 50 && y < 150) || (x > 300 && x < 400 && y > 100 && y < 200) ||
            (x > 200 && x < 250 && y > 300 && y < 350) || (x > 400 && x < 450 && y > 400 && y < 450)) {
          pixel = 1;
        }

        if ((x > 150 && x < 300 && y > 120 && y < 130)) {
          pixel = 1;
        }

        test_image.data[idx] = pixel;
      }
    }

    input_data_ = 1;  // Произвольное положительное число
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() < 3) {
      return false;
    }

    int num_components = output_data[2];
    return num_components > 0 && num_components <= 10;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(IvanovaPRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, IvanovaPMarkingComponentsOnBinaryImageSEQ>(
    PPC_SETTINGS_ivanova_p_marking_components_on_binary_image);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = IvanovaPRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, IvanovaPRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ivanova_p_marking_components_on_binary_image
