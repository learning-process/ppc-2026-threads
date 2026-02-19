#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "otcheskov_s_contrast_lin_stretch/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace otcheskov_s_contrast_lin_stretch {
namespace {
std::vector<uint8_t> CreateLowContrastImage(size_t size, uint8_t low = 100, uint8_t range = 50) {
  std::vector<uint8_t> image(size * size);
  for (size_t row = 0; row < size; ++row) {
    for (size_t col = 0; col < size; ++col) {
      uint8_t value = low + ((row + col) % range);
      image[(row * size) + col] = value;
    }
  }
  return image;
}

}  // namespace

class OtcheskovSContrastLinStretchPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kMatrixSize = 10000;
  InType input_img_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> low_dist(0, 200);
    int low = low_dist(gen);

    std::uniform_int_distribution<int> range_dist(5, 255 - low);
    int range = range_dist(gen);
    input_img_ = CreateLowContrastImage(kMatrixSize, low, range);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    return (*min_it == 0 && *max_it == 255) || (*min_it == *max_it);
  }

  InType GetTestInputData() final {
    return input_img_;
  }
};

TEST_P(OtcheskovSContrastLinStretchPerfTests, RunPerfTests) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, OtcheskovSContrastLinStretchSEQ>(PPC_SETTINGS_otcheskov_s_contrast_lin_stretch);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = OtcheskovSContrastLinStretchPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfTests, OtcheskovSContrastLinStretchPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace otcheskov_s_contrast_lin_stretch
