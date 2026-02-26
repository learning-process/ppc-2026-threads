#include <gtest/gtest.h>

#include <algorithm>
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

    int img_size = static_cast<int>(std::sqrt(size));
    if (img_size * img_size < size) {
      img_size++;
    }

    input_data_.width = img_size;
    input_data_.height = img_size;
    input_data_.data.resize(input_data_.width * input_data_.height);

    for (int i = 0; i < input_data_.width * input_data_.height; ++i) {
      input_data_.data[i] = (i % 101);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.width != input_data_.width || output_data.height != input_data_.height ||
        output_data.data.size() != input_data_.data.size()) {
      return false;
    }

    if (input_data_.width < 3 || input_data_.height < 3) {
      return true;
    }

    int cx = input_data_.width / 2;
    int cy = input_data_.height / 2;

    float expected = 0.0f;
    for (int ky = -1; ky <= 1; ++ky) {
      for (int kx = -1; kx <= 1; ++kx) {
        int px = cx + kx;
        int py = cy + ky;
        if (px < 0) {
          px = 0;
        }
        if (px >= input_data_.width) {
          px = input_data_.width - 1;
        }
        if (py < 0) {
          py = 0;
        }
        if (py >= input_data_.height) {
          py = input_data_.height - 1;
        }

        int kernel_idx = (ky + 1) * 3 + (kx + 1);
        expected += input_data_.data[py * input_data_.width + px] * kGaussKernel[kernel_idx];
      }
    }

    int actual = output_data.data[cy * output_data.width + cx];
    int expected_int = static_cast<int>(expected + 0.5f);

    return std::abs(actual - expected_int) <= 1;
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

INSTANTIATE_TEST_SUITE_P(GaussFilterTests, TerekhovDRunFuncTestsGauss, kGtestValues, kTestName);  // 1

}  // namespace

}  // namespace terekhov_d_seq_gauss_vert
