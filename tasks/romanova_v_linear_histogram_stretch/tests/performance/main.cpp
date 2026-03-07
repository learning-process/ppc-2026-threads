#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "romanova_v_linear_histogram_stretch/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

class RomanovaVRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const size_t kImageSize_ = 20000;
  InType input_data_;

  static std::vector<uint8_t> MakeImg(size_t width, size_t height, uint8_t low = 75, uint8_t range = 50) {
    std::vector<uint8_t> img(width * height);
    for (size_t i = 0; i < width; i++) {
      for (size_t j = 0; j < height; j++) {
        img[(i * height) + j] = low + ((j * 13 + i * 109) % range);
      }
    }
    return img;
  }

  void SetUp() override {
    input_data_ = MakeImg(kImageSize_, kImageSize_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    return (*min_it == 0 && *max_it == 255) || (*min_it == *max_it);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(RomanovaVRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RomanovaVLinHistogramStretchSEQ>(
    PPC_SETTINGS_romanova_v_linear_histogram_stretch);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RomanovaVRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RomanovaVRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace romanova_v_linear_histogram_stretch_threads
