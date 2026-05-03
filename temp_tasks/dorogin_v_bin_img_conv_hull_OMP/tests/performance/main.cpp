#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_OMP/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull_OMP/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorogin_v_bin_img_conv_hull_omp {

namespace {

BinaryImage MakeCheckerboard(int w, int h, int block) {
  BinaryImage img;
  img.width = w;
  img.height = h;
  img.data.assign(static_cast<std::size_t>(w) * static_cast<std::size_t>(h), 0);

  for (int yy = 0; yy < h; ++yy) {
    for (int xx = 0; xx < w; ++xx) {
      const bool on = ((xx / block) % 2 == 0) && ((yy / block) % 2 == 0);
      img.data[(static_cast<std::size_t>(yy) * static_cast<std::size_t>(w)) + static_cast<std::size_t>(xx)] =
          on ? 1U : 0U;
    }
  }

  return img;
}

}  // namespace

class DoroginVRunPerfTestsOMP : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    input_data_ = MakeCheckerboard(512, 512, 8);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty() || output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(DoroginVRunPerfTestsOMP, RunPerfOMP) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DoroginVBinImgConvHullOMP>(PPC_SETTINGS_dorogin_v_bin_img_conv_hull_OMP);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DoroginVRunPerfTestsOMP::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTests, DoroginVRunPerfTestsOMP, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull_omp
