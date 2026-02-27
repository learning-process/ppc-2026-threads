#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "badanov_a_select_edge_sobel/common/include/common.hpp"
#include "badanov_a_select_edge_sobel/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace badanov_a_select_edge_sobel {

class BadanovASelectEdgeSobelPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kWidth = 3840;
  static constexpr int kHeight = 2160;
  
  void SetUp() override {
    input_data_.resize(kWidth * kHeight);
    
    for (int row = 0; row < kHeight; ++row) {
      for (int col = 0; col < kWidth; ++col) {
        double val1 = std::sin(2.0 * std::numbers::pi * col / 32.0) * 
                     std::sin(2.0 * std::numbers::pi * row / 32.0);
        double val2 = std::sin(2.0 * std::numbers::pi * col / 16.0) * 
                     std::cos(2.0 * std::numbers::pi * row / 16.0);
        double val3 = std::sin(2.0 * std::numbers::pi * (col + row) / 24.0);
        
        double val = (val1 + val2 + val3) / 3.0;
        input_data_[row * kWidth + col] = static_cast<uint8_t>((val + 1.0) * 0.5 * 255.0);
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != static_cast<size_t>(kWidth * kHeight)) {
      return false;
    }
    
    for (int col = 0; col < kWidth; ++col) {
      if (output_data[col] != 0) return false;
      if (output_data[((kHeight - 1) * kWidth) + col] != 0) return false;
    }
    
    for (int row = 0; row < kHeight; ++row) {
      if (output_data[row * kWidth] != 0) return false;
      if (output_data[(row * kWidth) + (kWidth - 1)] != 0) return false;
    }
    
    bool has_edges = false;
    for (int row = 1; row < kHeight - 1 && !has_edges; ++row) {
      for (int col = 1; col < kWidth - 1 && !has_edges; ++col) {
        if (output_data[(row * kWidth) + col] > 0) {
          has_edges = true;
        }
      }
    }
    
    return has_edges;
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
    ppc::util::MakeAllPerfTasks<InType, BadanovASelectEdgeSobel>(PPC_SETTINGS_badanov_a_select_edge_sobel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BadanovASelectEdgeSobelPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BadanovASelectEdgeSobelPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace badanov_a_select_edge_sobel