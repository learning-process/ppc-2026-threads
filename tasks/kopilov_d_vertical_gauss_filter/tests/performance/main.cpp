#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class GaussianFilterPerformanceTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kImageWidth = 8192;
  static constexpr int kImageHeight = 8192;
  InType inputImage_{};

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    inputImage_.width = kImageWidth;
    inputImage_.height = kImageHeight;

    inputImage_.data.resize(static_cast<size_t>(kImageWidth) * static_cast<size_t>(kImageHeight));
    for (auto &val : inputImage_.data) {
      val = static_cast<std::uint8_t>(dist(gen));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.width == inputImage_.width && output_data.height == inputImage_.height &&
           output_data.data.size() == inputImage_.data.size();
  }

  InType GetTestInputData() final {
    return inputImage_;
  }
};

TEST_P(GaussianFilterPerformanceTests, PerformanceTest) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KopilovDVerticalGaussFilterSEQ, KopilovDVerticalGaussFilterOMP>(
        PPC_SETTINGS_kopilov_d_vertical_gauss_filter);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GaussianFilterPerformanceTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GaussianFilterPerformanceTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kopilov_d_vertical_gauss_filter
