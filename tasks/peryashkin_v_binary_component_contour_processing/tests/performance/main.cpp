// tasks/peryashkin_v_binary_component_contour_processing/tests/performance/main.cpp
#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <vector>

#include "peryashkin_v_binary_component_contour_processing/common/include/common.hpp"
#include "peryashkin_v_binary_component_contour_processing/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace peryashkin_v_binary_component_contour_processing {

namespace {

BinaryImage MakePattern(int w, int h, int step) {
  BinaryImage img;
  img.width = w;
  img.height = h;
  img.data.assign(static_cast<std::size_t>(w) * static_cast<std::size_t>(h), 0);

  for (int yy = 0; yy < h; ++yy) {
    for (int xx = 0; xx < w; ++xx) {
      const bool on = ((xx / step) % 2 == 0) && ((yy / step) % 2 == 0);
      img.data[(static_cast<std::size_t>(yy) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(xx)] =
          on ? 1 : 0;
    }
  }
  return img;
}

}  // namespace

class PeryashkinVRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    input_data_ = MakePattern(512, 512, 8);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    (void)output_data;
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PeryashkinVRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, PeryashkinVBinaryComponentContourProcessingSEQ>(
    PPC_SETTINGS_peryashkin_v_binary_component_contour_processing);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PeryashkinVRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTests, PeryashkinVRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace peryashkin_v_binary_component_contour_processing
