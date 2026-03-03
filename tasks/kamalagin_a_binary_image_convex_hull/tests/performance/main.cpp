#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace kamalagin_a_binary_image_convex_hull {

namespace {

BinaryImage MakeSyntheticPerfImage() {
  const int rows = 100;
  const int cols = 100;
  BinaryImage img;
  img.rows = rows;
  img.cols = cols;
  img.data.resize(static_cast<size_t>(rows * cols), 0);
  for (int r = 10; r < 30; ++r) {
    for (int c = 10; c < 40; ++c) {
      img.data[static_cast<size_t>(r * cols + c)] = 1;
    }
  }
  for (int r = 50; r < 70; ++r) {
    for (int c = 40; c < 70; ++c) {
      img.data[static_cast<size_t>(r * cols + c)] = 1;
    }
  }
  for (int r = 80; r < 95; ++r) {
    for (int c = 5; c < 35; ++c) {
      img.data[static_cast<size_t>(r * cols + c)] = 1;
    }
  }
  for (int i = 0; i < rows; ++i) {
    int c = i % cols;
    img.data[static_cast<size_t>(i * cols + c)] = 1;
  }
  return img;
}

}  // namespace

class KamalaginRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    input_data_ = MakeSyntheticPerfImage();
  }

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::steady_clock::now();
    perf_attrs.current_timer = [t0] {
      auto now = std::chrono::steady_clock::now();
      auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - t0).count();
      return static_cast<double>(ns) * 1e-9;
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  InType input_data_{};
};

TEST_P(KamalaginRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KamalaginABinaryImageConvexHullSEQ>(
    PPC_SETTINGS_kamalagin_a_binary_image_convex_hull);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = KamalaginRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KamalaginRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamalagin_a_binary_image_convex_hull
