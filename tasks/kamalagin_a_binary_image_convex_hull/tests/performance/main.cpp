#include <gtest/gtest.h>

#include <chrono>
#include <cstddef>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/all/include/ops_all.hpp"
#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"
#include "kamalagin_a_binary_image_convex_hull/omp/include/ops_omp.hpp"
#include "kamalagin_a_binary_image_convex_hull/seq/include/ops_seq.hpp"
#include "kamalagin_a_binary_image_convex_hull/stl/include/ops_stl.hpp"
#include "kamalagin_a_binary_image_convex_hull/tbb/include/ops_tbb.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace kamalagin_a_binary_image_convex_hull {

namespace {

BinaryImage MakeSyntheticPerfImage() {
  constexpr int kRows = 5000;
  constexpr int kCols = 5000;
  BinaryImage img;
  img.rows = kRows;
  img.cols = kCols;
  const size_t total = static_cast<size_t>(kRows) * static_cast<size_t>(kCols);
  img.data.assign(total, 0);
  for (int row = 0; row < kRows; ++row) {
    const auto urow = static_cast<size_t>(row);
    const int col1 = (row * 13) % kCols;
    const int col2 = (row * 29 + 7) % kCols;
    img.data[(urow * static_cast<size_t>(kCols)) + static_cast<size_t>(col1)] = 1;
    img.data[(urow * static_cast<size_t>(kCols)) + static_cast<size_t>(col2)] = 1;
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

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KamalaginABinaryImageConvexHullSEQ, KamalaginABinaryImageConvexHullOMP,
                                KamalaginABinaryImageConvexHullTBB, KamalaginABinaryImageConvexHullSTL,
                                KamalaginABinaryImageConvexHullALL>(PPC_SETTINGS_kamalagin_a_binary_image_convex_hull);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = KamalaginRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KamalaginRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamalagin_a_binary_image_convex_hull
