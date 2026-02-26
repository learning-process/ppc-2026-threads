#include <gtest/gtest.h>

#include <cstdint>
#include <random>
#include <vector>

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "iskhakov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace iskhakov_d_vertical_gauss_filter {

class IskhakovDVerticalGaussFilterPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kWidth = 4096;
  const int kHeight = 4096;
  InType input_data_{};

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    input_data_.width = kWidth;
    input_data_.height = kHeight;
    input_data_.data.resize(kWidth * kHeight);
    for (auto &val : input_data_.data) {
      val = static_cast<std::uint8_t>(dist(gen));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.width == input_data_.width && output_data.height == input_data_.height &&
           output_data.data.size() == input_data_.data.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(IskhakovDVerticalGaussFilterPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, IskhakovDVerticalGaussFilterSEQ>(PPC_SETTINGS_iskhakov_d_vertical_gauss_filter);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = IskhakovDVerticalGaussFilterPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, IskhakovDVerticalGaussFilterPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace iskhakov_d_vertical_gauss_filter
