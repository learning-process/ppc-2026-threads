#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <random>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_TBB/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull_TBB/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorogin_v_bin_img_conv_hull_tbb {

class DoroginVImgConvHullTBBPerformanceTest : public ppc::util::BaseRunPerfTests<InputType, OutputType> {
  static constexpr int kImageSize = 600;

  void SetUp() override {
    input_.rows = kImageSize;
    input_.cols = kImageSize;
    input_.pixels.assign(static_cast<size_t>(kImageSize) * static_cast<size_t>(kImageSize), 0);

    std::mt19937 rng(static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::uniform_int_distribution<int> dist(0, kImageSize - 1);

    for (int i = 0; i < 5; ++i) {
      const int x0 = dist(rng) % (kImageSize - 35);
      const int y0 = dist(rng) % (kImageSize - 35);
      const int x1 = x0 + 20 + (dist(rng) % 15);
      const int y1 = y0 + 20 + (dist(rng) % 15);
      for (int row_idx = y0; row_idx <= y1; ++row_idx) {
        for (int col_idx = x0; col_idx <= x1; ++col_idx) {
          input_
              .pixels[(static_cast<size_t>(row_idx) * static_cast<size_t>(kImageSize)) + static_cast<size_t>(col_idx)] =
              255;
        }
      }
    }

    for (int i = 0; i < kImageSize; i += 19) {
      input_.pixels[(static_cast<size_t>(i) * static_cast<size_t>(kImageSize)) + static_cast<size_t>(i)] = 255;
    }
  }

  bool CheckTestOutputData(OutputType &output) final {
    if (output.empty()) {
      return false;
    }
    return std::ranges::all_of(output, [](const auto &hull) { return !hull.empty(); });
  }

  InputType GetTestInputData() final {
    return input_;
  }

  InputType input_;
};

TEST_P(DoroginVImgConvHullTBBPerformanceTest, DoroginVPerfBinaryImageConvexHullTBB) {
  ExecuteTest(GetParam());
}

namespace {
const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InputType, DoroginVImgConvHullTBB>(PPC_SETTINGS_dorogin_v_bin_img_conv_hull_TBB);
const auto kValues = ppc::util::TupleToGTestValues(kPerfTasks);
INSTANTIATE_TEST_SUITE_P(DoroginVBinImgConvHullTBBPerformance, DoroginVImgConvHullTBBPerformanceTest, kValues,
                         DoroginVImgConvHullTBBPerformanceTest::CustomPerfTestName);
}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull_tbb
