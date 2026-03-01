#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "egorova_l_binary_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace egorova_l_binary_convex_hull {

class EgorovaLPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kWidth = 1000;
  static constexpr int kHeight = 1000;
  InType input_data_;
  int expected_components_{};

  void SetUp() override {
    input_data_.width = kWidth;
    input_data_.height = kHeight;
    const size_t image_size = static_cast<size_t>(kWidth) * static_cast<size_t>(kHeight);
    input_data_.data.assign(image_size, 0);

    // Draw first rectangle
    for (int row = 100; row < 400; ++row) {
      for (int col = 100; col < 400; ++col) {
        input_data_.data[(static_cast<size_t>(row) * static_cast<size_t>(kWidth)) + static_cast<size_t>(col)] = 255;
      }
    }

    // Draw second rectangle
    for (int row = 500; row < 700; ++row) {
      for (int col = 500; col < 700; ++col) {
        input_data_.data[(static_cast<size_t>(row) * static_cast<size_t>(kWidth)) + static_cast<size_t>(col)] = 255;
      }
    }

    // Draw third rectangle
    for (int row = 800; row < 850; ++row) {
      for (int col = 800; col < 850; ++col) {
        input_data_.data[(static_cast<size_t>(row) * static_cast<size_t>(kWidth)) + static_cast<size_t>(col)] = 255;
      }
    }

    expected_components_ = 3;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != static_cast<size_t>(expected_components_)) {
      return false;
    }

    return std::ranges::all_of(output_data, [](const std::vector<Point> &hull) { return hull.size() >= 3; });
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(EgorovaLPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BinaryConvexHullSEQ>(PPC_SETTINGS_egorova_l_binary_convex_hull);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = EgorovaLPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, EgorovaLPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace egorova_l_binary_convex_hull
