#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <vector>

#include "terekhov_d_seq_gauss_vert/common/include/common.hpp"
#include "terekhov_d_seq_gauss_vert/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace terekhov_d_seq_gauss_vert {

class TerekhovDGaussVertSEQPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const size_t total_pixels = 10000000; // 10M пикселей
    int img_size = static_cast<int>(std::sqrt(total_pixels));
    
    input_data_.width = img_size;
    input_data_.height = img_size;
    input_data_.data.resize(input_data_.width * input_data_.height);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);
    
    for (auto &pixel : input_data_.data) {
      pixel = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.width == input_data_.width && 
           output_data.height == input_data_.height &&
           output_data.data.size() == input_data_.data.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(TerekhovDGaussVertSEQPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TerekhovDGaussVertSEQ>(PPC_SETTINGS_terekhov_d_seq_gauss_vert);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TerekhovDGaussVertSEQPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TerekhovDGaussVertSEQPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace terekhov_d_seq_gauss_vert