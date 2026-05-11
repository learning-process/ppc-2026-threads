#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "dorogin_v_bin_img_conv_hull_all/all/include/ops_all.hpp"
#include "dorogin_v_bin_img_conv_hull_all/common/include/common.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorogin_v_bin_img_conv_hull_all {

class DoroginVRunPerfTestsALL : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    constexpr int kWidth = 256;
    constexpr int kHeight = 256;
    std::vector<std::uint8_t> data(static_cast<std::size_t>(kWidth) * static_cast<std::size_t>(kHeight), 0U);
    for (int row = 24; row < 120; ++row) {
      for (int col = 18; col < 125; ++col) {
        data[(static_cast<std::size_t>(row) * static_cast<std::size_t>(kWidth)) + static_cast<std::size_t>(col)] = 1U;
      }
    }
    for (int row = 130; row < 230; ++row) {
      for (int col = 145; col < 245; ++col) {
        data[(static_cast<std::size_t>(row) * static_cast<std::size_t>(kWidth)) + static_cast<std::size_t>(col)] = 1U;
      }
    }
    input_data_ = InType{.width = kWidth, .height = kHeight, .data = std::move(data)};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(DoroginVRunPerfTestsALL, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DoroginVConvHullAll>(PPC_SETTINGS_dorogin_v_bin_img_conv_hull_all);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DoroginVRunPerfTestsALL::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DoroginVRunPerfTestsALL, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dorogin_v_bin_img_conv_hull_all
