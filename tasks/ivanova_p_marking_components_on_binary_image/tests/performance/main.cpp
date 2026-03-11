#include <gtest/gtest.h>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"
#include "ivanova_p_marking_components_on_binary_image/data/image_generator.hpp"
#include "ivanova_p_marking_components_on_binary_image/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace ivanova_p_marking_components_on_binary_image {

class IvanovaPRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 512;
  InType input_data_ = 0;

  void SetUp() override {
    test_image = CreateTestImage(kSize_, kSize_, 8);  // Например, тест 8 с множеством компонент

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
