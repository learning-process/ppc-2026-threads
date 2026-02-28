#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <vector>

#include "marin_l_mark_components_seq/common/include/common.hpp"
#include "marin_l_mark_components_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace marin_l_mark_components_seq {

namespace {

Image MakeRandomBinaryImage(int height, int width, double fill_probability) {
  Image image(static_cast<size_t>(height), std::vector<int>(static_cast<size_t>(width), 0));

  std::random_device random_seed;
  std::mt19937 generator(random_seed());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for (int row_idx = 0; row_idx < height; ++row_idx) {
    for (int col_idx = 0; col_idx < width; ++col_idx) {
      image[static_cast<size_t>(row_idx)][static_cast<size_t>(col_idx)] =
          (distribution(generator) < fill_probability) ? 1 : 0;
    }
  }
  return image;
}

}  // namespace

class MarinLRunPerfTestComponents : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_width_pixels = 1024;
  const int k_height_pixels = 1024;
  InType input_data{};

  void SetUp() override {
    input_data.binary = MakeRandomBinaryImage(k_height_pixels, k_width_pixels, 0.3);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &output_labels = output_data.labels;

    if (output_labels.size() != input_data.binary.size()) {
      return false;
    }
    if (!output_labels.empty() && output_labels[0].size() != input_data.binary[0].size()) {
      return false;
    }

    const int height_pixels = static_cast<int>(output_labels.size());
    const int width_pixels = height_pixels != 0 ? static_cast<int>(output_labels[0].size()) : 0;

    for (int row_idx = 0; row_idx < height_pixels; ++row_idx) {
      for (int col_idx = 0; col_idx < width_pixels; ++col_idx) {
        if (input_data.binary[row_idx][col_idx] == 0 && output_labels[row_idx][col_idx] != 0) {
          return false;
        }
        if (input_data.binary[row_idx][col_idx] == 1 && output_labels[row_idx][col_idx] <= 0) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(MarinLRunPerfTestComponents, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MarinLMarkComponentsSEQ>(PPC_SETTINGS_marin_l_mark_components_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MarinLRunPerfTestComponents::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ComponentLabelingPerf, MarinLRunPerfTestComponents, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace marin_l_mark_components_seq
