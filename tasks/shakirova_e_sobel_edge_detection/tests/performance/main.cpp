#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "shakirova_e_sobel_edge_detection/common/include/common.hpp"
#include "shakirova_e_sobel_edge_detection/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shakirova_e_sobel_edge_detection {

class ShakirovaESobelEdgeDetectionPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kWidth = 3840;
  static constexpr int kHeight = 2160;
  InType input_data_;

  void SetUp() override {
    input_data_ = ImgContainer(kWidth, kHeight);

    for (int row = 0; row < kHeight; ++row) {
      for (int col = 0; col < kWidth; ++col) {
        const double val = std::sin(2.0 * std::numbers::pi * col / 8.0) * std::sin(2.0 * std::numbers::pi * row / 8.0);
        input_data_.pixels[(row * kWidth) + col] = static_cast<int>((val + 1.0) * 0.5 * 255.0);
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (static_cast<int>(output_data.size()) != kWidth * kHeight) {
      return false;
    }

    for (int col = 0; col < kWidth; ++col) {
      if (output_data[col] != 0) {
        return false;
      }
      if (output_data[((kHeight - 1) * kWidth) + col] != 0) {
        return false;
      }
    }

    int edge_count = 0;
    for (int row = 1; row < kHeight - 1; ++row) {
      for (int col = 1; col < kWidth - 1; ++col) {
        if (output_data[(row * kWidth) + col] > 0) {
          ++edge_count;
        }
      }
    }
    const int inner = (kWidth - 2) * (kHeight - 2);
    return edge_count >= inner / 2;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ShakirovaESobelEdgeDetectionPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShakirovaESobelEdgeDetectionSEQ>(PPC_SETTINGS_shakirova_e_sobel_edge_detection);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShakirovaESobelEdgeDetectionPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShakirovaESobelEdgeDetectionPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shakirova_e_sobel_edge_detection
