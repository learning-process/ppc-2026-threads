#include <gtest/gtest.h>

#include <chrono>
#include <cstddef>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull/seq/include/ops_seq.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace kamalagin_a_binary_image_convex_hull {

namespace {

BinaryImage MakeSyntheticPerfImage() {
  const int rows = 100;
  const int cols = 100;
  BinaryImage img;
  img.rows = rows;
  img.cols = cols;
  img.data.resize(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0);
  for (int row = 10; row < 30; ++row) {
    for (int col = 10; col < 40; ++col) {
      img.data[(static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col)] = 1;
    }
  }
  for (int row = 50; row < 70; ++row) {
    for (int col = 40; col < 70; ++col) {
      img.data[(static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col)] = 1;
    }
  }
  for (int row = 80; row < 95; ++row) {
    for (int col = 5; col < 35; ++col) {
      img.data[(static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col)] = 1;
    }
  }
  for (int i = 0; i < rows; ++i) {
    int col = i % cols;
    img.data[(static_cast<size_t>(i) * static_cast<size_t>(cols)) + static_cast<size_t>(col)] = 1;
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
