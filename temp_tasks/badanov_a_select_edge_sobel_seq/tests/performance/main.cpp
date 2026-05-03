#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "badanov_a_select_edge_sobel_seq/common/include/common.hpp"
#include "badanov_a_select_edge_sobel_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace badanov_a_select_edge_sobel_seq {

class BadanovASelectEdgeSobelPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kWidth = 3840;
  static constexpr int kHeight = 2160;

  void SetUp() override {
    const size_t total_pixels = static_cast<size_t>(kWidth) * static_cast<size_t>(kHeight);
    input_data_.resize(total_pixels);

    for (int row = 0; row < kHeight; ++row) {
      for (int col = 0; col < kWidth; ++col) {
        const double value = static_cast<double>(col % 256) / 255.0;
        const size_t index = (static_cast<size_t>(row) * static_cast<size_t>(kWidth)) + static_cast<size_t>(col);
        input_data_[index] = static_cast<uint8_t>(value * 255.0);
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const size_t expected_size = static_cast<size_t>(kWidth) * static_cast<size_t>(kHeight);
    return output_data.size() == expected_size;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(BadanovASelectEdgeSobelPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BadanovASelectEdgeSobelSEQ>(PPC_SETTINGS_badanov_a_select_edge_sobel_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BadanovASelectEdgeSobelPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BadanovASelectEdgeSobelPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace badanov_a_select_edge_sobel_seq
