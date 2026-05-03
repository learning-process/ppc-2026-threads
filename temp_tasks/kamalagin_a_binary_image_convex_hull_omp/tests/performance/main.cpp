#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull_omp/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace kamalagin_a_binary_image_convex_hull_omp {

namespace {

BinaryImage MakePerfImage() {
  constexpr int kSize = 100;
  std::vector<uint8_t> data(static_cast<size_t>(kSize) * static_cast<size_t>(kSize), 0);
  for (int row = 10; row < 40; ++row) {
    for (int col = 10; col < 40; ++col) {
      data[(static_cast<size_t>(row) * static_cast<size_t>(kSize)) + static_cast<size_t>(col)] = 1;
    }
  }
  for (int row = 50; row < 90; ++row) {
    for (int col = 50; col < 90; ++col) {
      data[(static_cast<size_t>(row) * static_cast<size_t>(kSize)) + static_cast<size_t>(col)] = 1;
    }
  }
  for (int i = 0; i < kSize; ++i) {
    data[(static_cast<size_t>(i) * static_cast<size_t>(kSize)) + static_cast<size_t>(i)] = 1;
  }
  return BinaryImage{.rows = kSize, .cols = kSize, .data = std::move(data)};
}

}  // namespace

class KamalaginABinaryImageConvexHullPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    input_data_ = MakePerfImage();
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KamalaginABinaryImageConvexHullPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KamalaginABinaryImageConvexHullOMP>(
    PPC_SETTINGS_kamalagin_a_binary_image_convex_hull_omp);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KamalaginABinaryImageConvexHullPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KamalaginABinaryImageConvexHullPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamalagin_a_binary_image_convex_hull_omp
