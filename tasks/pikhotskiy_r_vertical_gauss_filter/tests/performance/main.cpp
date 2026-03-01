#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "pikhotskiy_r_vertical_gauss_filter/common/include/common.hpp"
#include "pikhotskiy_r_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace pikhotskiy_r_vertical_gauss_filter {

class PikhotskiyRVerticalGaussFilterPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kImgWidth = 4096;
  static constexpr int kImgHeight = 4096;
  InType input_img_{};

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    input_img_.width = kImgWidth;
    input_img_.height = kImgHeight;

    input_img_.data.resize(static_cast<size_t>(kImgWidth) * static_cast<size_t>(kImgHeight));
    for (auto &val : input_img_.data) {
      val = static_cast<std::uint8_t>(dist(gen));
    }
  }

  bool CheckTestOutputData(OutType &output) final {
    return output.width == input_img_.width && output.height == input_img_.height &&
           output.data.size() == input_img_.data.size();
  }

  InType GetTestInputData() final {
    return input_img_;
  }
};

TEST_P(PikhotskiyRVerticalGaussFilterPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PikhotskiyRVerticalGaussFilterSEQ>(PPC_SETTINGS_pikhotskiy_r_vertical_gauss_filter);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PikhotskiyRVerticalGaussFilterPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfRunTests, PikhotskiyRVerticalGaussFilterPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pikhotskiy_r_vertical_gauss_filter