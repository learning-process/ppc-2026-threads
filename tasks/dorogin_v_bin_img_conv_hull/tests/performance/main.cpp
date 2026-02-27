#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "dorogin_v_bin_img_conv_hull/common/include/common.hpp"
#include "dorogin_v_bin_img_conv_hull/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorogin_v_bin_img_conv_hull {

class DoroginVBinImgConvHullPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kSize = 600;

 protected:
  void SetUp() override {
    input_.width = static_cast<int>(kSize);
    input_.height = static_cast<int>(kSize);
    input_.pixels.assign(kSize * kSize, 0);

    // Диагональные компоненты + псевдослучайные точки
    for (size_t i = 0; i < kSize; ++i) {
      size_t idx1 = i * kSize + i;
      size_t idx2 = i * kSize + (kSize - i - 1);

      input_.pixels[idx1] = 255;
      input_.pixels[idx2] = 255;

      if (i % 17 == 0) {
        size_t extra = i * kSize + ((i * 37 + 11) % kSize);
        input_.pixels[extra] = 255;
      }
    }
  }

  bool CheckTestOutputData(OutType &out) override {
    return !out.convex_hulls.empty();
  }

  InType GetTestInputData() override {
    return input_;
  }

 private:
  InType input_;
};

TEST_P(DoroginVBinImgConvHullPerfTests, RunPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, doroginVBinImgConvHullSeq>(PPC_SETTINGS_dorogin_v_bin_img_conv_hull);

const auto kValues = ppc::util::TupleToGTestValues(kPerfTasks);

INSTANTIATE_TEST_SUITE_P(Perf, DoroginVBinImgConvHullPerfTests, kValues,
                         DoroginVBinImgConvHullPerfTests::CustomPerfTestName);

}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull
