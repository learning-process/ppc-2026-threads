#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "terekhov_d_seq_gauss_vert/common/include/common.hpp"
#include "terekhov_d_seq_gauss_vert/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace terekhov_d_seq_gauss_vert {

class TerekhovDRunFuncTestsGauss : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

 protected:
  void SetUp() override {
    TestType size = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    int img_size = static_cast<int>(std::sqrt(static_cast<double>(size)));
    if (img_size * img_size < size) {
      ++img_size;
    }

    input_data_.width = img_size;
    input_data_.height = img_size;
    input_data_.data.resize(static_cast<size_t>(input_data_.width) * static_cast<size_t>(input_data_.height));

    for (size_t i = 0; i < input_data_.data.size(); ++i) {
      input_data_.data[i] = static_cast<int>(i % 101);
    }
  }

bool CheckTestOutputData(OutType &output_data) final {
  if (!ValidateOutputSize(output_data)) {
    return false;
  }
  
  if (input_data_.width < 3 || input_data_.height < 3) {
    return true;
  }
  
  return ValidateCenterPixel(output_data);
}

private:
  bool ValidateOutputSize(const OutType &output_data) const {
    if (output_data.width != input_data_.width || 
        output_data.height != input_data_.height) {
      return false;
    }
    
    if (output_data.data.size() != input_data_.data.size()) {
      return false;
    }
    
    return true;
  }
  
  bool ValidateCenterPixel(const OutType &output_data) const {
    int cx = input_data_.width / 2;
    int cy = input_data_.height / 2;
    
    float expected = ComputeExpectedValue(cx, cy);
    int actual = GetActualValue(output_data, cx, cy);
    int expected_int = static_cast<int>(std::lround(expected));
    
    return std::abs(actual - expected_int) <= 1;
  }
  
  float ComputeExpectedValue(int cx, int cy) const {
    float sum = 0.0F;
    
    for (int ky = -1; ky <= 1; ++ky) {
      for (int kx = -1; kx <= 1; ++kx) {
        int px = cx + kx;
        int py = cy + ky;
        
        px = ClampCoordinate(px, 0, input_data_.width - 1);
        py = ClampCoordinate(py, 0, input_data_.height - 1);
        
        int kernel_idx = ((ky + 1) * 3) + (kx + 1);
        size_t data_idx = (static_cast<size_t>(py) * static_cast<size_t>(input_data_.width)) + static_cast<size_t>(px);
        
        sum += static_cast<float>(input_data_.data[data_idx]) * kGaussKernel[static_cast<size_t>(kernel_idx)];
      }
    }
    
    return sum;
  }
  
  int GetActualValue(const OutType &output_data, int cx, int cy) const {
    size_t out_idx = (static_cast<size_t>(cy) * static_cast<size_t>(output_data.width)) + static_cast<size_t>(cx);
    return output_data.data[out_idx];
  }
  
  int ClampCoordinate(int value, int min_val, int max_val) const {
    if (value < min_val) return min_val;
    if (value > max_val) return max_val;
    return value;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(TerekhovDRunFuncTestsGauss, GaussFilter) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {16, 256, 1024};

const auto kTestTasksList =
    ppc::util::AddFuncTask<TerekhovDGaussVertSEQ, InType>(kTestParam, PPC_SETTINGS_terekhov_d_seq_gauss_vert);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = TerekhovDRunFuncTestsGauss::PrintFuncTestName<TerekhovDRunFuncTestsGauss>;

INSTANTIATE_TEST_SUITE_P(GaussFilterTests, TerekhovDRunFuncTestsGauss, kGtestValues, kTestName);

}  // namespace

}  // namespace terekhov_d_seq_gauss_vert
