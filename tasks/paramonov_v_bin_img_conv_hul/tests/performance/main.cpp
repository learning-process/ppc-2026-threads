#include <gtest/gtest.h>

#include <chrono>
#include <random>
#include <vector>

#include "paramonov_v_bin_img_conv_hul/common/include/common.hpp"
#include "paramonov_v_bin_img_conv_hul/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace paramonov_v_bin_img_conv_hul {

class ConvexHullPerformanceTest : public ppc::util::BaseRunPerfTests<InputType, OutputType> {
  static constexpr int kImageSize = 600;

  void SetUp() override {
    input_image_.rows = kImageSize;
    input_image_.cols = kImageSize;
    input_image_.pixels.assign(static_cast<size_t>(kImageSize) * kImageSize, 0);

    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> dist(0, kImageSize - 1);

    // Создаем 5 случайных прямоугольников
    for (int comp = 0; comp < 5; ++comp) {
      int x1 = dist(rng) % (kImageSize - 30);
      int y1 = dist(rng) % (kImageSize - 30);
      int x2 = x1 + 20 + dist(rng) % 20;
      int y2 = y1 + 20 + dist(rng) % 20;

      for (int r = y1; r <= y2; ++r) {
        for (int c = x1; c <= x2; ++c) {
          if (r >= 0 && r < kImageSize && c >= 0 && c < kImageSize) {
            size_t idx = static_cast<size_t>(r) * kImageSize + c;
            input_image_.pixels[idx] = 255;
          }
        }
      }
    }

    // Добавляем несколько диагональных линий
    for (int i = 0; i < kImageSize; i += 17) {
      int r = i;
      int c = i;
      if (r < kImageSize && c < kImageSize) {
        size_t idx = static_cast<size_t>(r) * kImageSize + c;
        input_image_.pixels[idx] = 255;
      }
    }
  }

  bool CheckTestOutputData(OutputType &output) override {
    if (output.empty()) {
      return false;
    }
    return std::all_of(output.begin(), output.end(), [](const auto &hull) { return !hull.empty(); });
  }

  InputType GetTestInputData() override {
    return input_image_;
  }

  InputType input_image_;
};

TEST_P(ConvexHullPerformanceTest, RunPerformanceTest) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerformanceTasks =
    ppc::util::MakeAllPerfTasks<InputType, ConvexHullSequential>(PPC_SETTINGS_paramonov_v_bin_img_conv_hul);

const auto kTestValues = ppc::util::TupleToGTestValues(kPerformanceTasks);

INSTANTIATE_TEST_SUITE_P(ParamonovPerfTests, ConvexHullPerformanceTest, kTestValues,
                         ConvexHullPerformanceTest::CustomPerfTestName);

}  // namespace

}  // namespace paramonov_v_bin_img_conv_hul
