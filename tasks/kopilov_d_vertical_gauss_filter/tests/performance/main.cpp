#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class KopilovDVerticalGaussFilterPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kWidth = 8192;
  static constexpr int kHeight = 8192;
  InType input_data_{};

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    input_data_.width = kWidth;
    input_data_.height = kHeight;

    input_data_.data.resize(static_cast<size_t>(kWidth) * static_cast<size_t>(kHeight));
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

TEST_P(KopilovDVerticalGaussFilterPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KopilovDVerticalGaussFilterSEQ>(PPC_SETTINGS_kopilov_d_vertical_gauss_filter);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KopilovDVerticalGaussFilterPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KopilovDVerticalGaussFilterPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kopilov_d_vertical_gauss_filter