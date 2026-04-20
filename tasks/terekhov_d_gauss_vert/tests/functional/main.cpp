#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "terekhov_d_gauss_vert/common/include/common.hpp"
#include "terekhov_d_gauss_vert/omp/include/ops_omp.hpp"
#include "terekhov_d_gauss_vert/seq/include/ops_seq.hpp"
#include "terekhov_d_gauss_vert/stl/include/ops_stl.hpp"
#include "terekhov_d_gauss_vert/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace terekhov_d_gauss_vert {

class TerekhovDGaussVertFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

  InType GetTestInputData() final {
    return input_data_;
  }

  [[nodiscard]] bool ValidateOutputSize(const OutType &output_data) const {
    return output_data.width == input_data_.width && output_data.height == input_data_.height &&
           output_data.data.size() == input_data_.data.size();
  }

  [[nodiscard]] bool ValidateCenterPixel(const OutType &output_data) const {
    int cx = input_data_.width / 2;
    int cy = input_data_.height / 2;
    float expected = ComputeExpectedValue(cx, cy);
    int actual = GetActualValue(output_data, cx, cy);
    return std::abs(actual - static_cast<int>(std::lround(expected))) <= 1;
  }

  [[nodiscard]] float ComputeExpectedValue(int cx, int cy) const {
    float sum = 0.0F;
    for (int ky = -1; ky <= 1; ++ky) {
      for (int kx = -1; kx <= 1; ++kx) {
        int px = ClampCoordinate(cx + kx, 0, input_data_.width - 1);
        int py = ClampCoordinate(cy + ky, 0, input_data_.height - 1);
        int kernel_idx = ((ky + 1) * 3) + (kx + 1);
        size_t data_idx = (static_cast<size_t>(py) * static_cast<size_t>(input_data_.width)) + static_cast<size_t>(px);
        sum += static_cast<float>(input_data_.data[data_idx]) * kGaussKernel[static_cast<size_t>(kernel_idx)];
      }
    }
    return sum;
  }

  [[nodiscard]] static int GetActualValue(const OutType &output_data, int cx, int cy) {
    size_t out_idx = (static_cast<size_t>(cy) * static_cast<size_t>(output_data.width)) + static_cast<size_t>(cx);
    return output_data.data[out_idx];
  }

  [[nodiscard]] static int ClampCoordinate(int value, int min_val, int max_val) {
    if (value < min_val) {
      return min_val;
    }
    if (value > max_val) {
      return max_val;
    }
    return value;
  }

 private:
  InType input_data_;
};

using TerekhovDGaussVertFuncTestsSEQ = TerekhovDGaussVertFuncTests;
using TerekhovDGaussVertFuncTestsOMP = TerekhovDGaussVertFuncTests;
using TerekhovDGaussVertFuncTestsTBB = TerekhovDGaussVertFuncTests;
using TerekhovDGaussVertFuncTestsSTL = TerekhovDGaussVertFuncTests;

namespace {

TEST_P(TerekhovDGaussVertFuncTestsSEQ, GaussFilterSEQ) {
  ExecuteTest(GetParam());
}
TEST_P(TerekhovDGaussVertFuncTestsOMP, GaussFilterOMP) {
  ExecuteTest(GetParam());
}
TEST_P(TerekhovDGaussVertFuncTestsTBB, GaussFilterTBB) {
  ExecuteTest(GetParam());
}
TEST_P(TerekhovDGaussVertFuncTestsSTL, GaussFilterSTL) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {16, 256, 1024};

const auto kSeqTasksList =
    ppc::util::AddFuncTask<TerekhovDGaussVertSEQ, InType>(kTestParam, PPC_SETTINGS_terekhov_d_gauss_vert);
const auto kOmpTasksList =
    ppc::util::AddFuncTask<TerekhovDGaussVertOMP, InType>(kTestParam, PPC_SETTINGS_terekhov_d_gauss_vert);
const auto kTbbTasksList =
    ppc::util::AddFuncTask<TerekhovDGaussVertTBB, InType>(kTestParam, PPC_SETTINGS_terekhov_d_gauss_vert);
const auto kStlTasksList =
    ppc::util::AddFuncTask<TerekhovDGaussVertSTL, InType>(kTestParam, PPC_SETTINGS_terekhov_d_gauss_vert);

const auto kSeqValues = ppc::util::ExpandToValues(kSeqTasksList);
const auto kOmpValues = ppc::util::ExpandToValues(kOmpTasksList);
const auto kTbbValues = ppc::util::ExpandToValues(kTbbTasksList);
const auto kStlValues = ppc::util::ExpandToValues(kStlTasksList);

const auto kNameFnSEQ = TerekhovDGaussVertFuncTestsSEQ::PrintFuncTestName<TerekhovDGaussVertFuncTestsSEQ>;
const auto kNameFnOMP = TerekhovDGaussVertFuncTestsOMP::PrintFuncTestName<TerekhovDGaussVertFuncTestsOMP>;
const auto kNameFnTBB = TerekhovDGaussVertFuncTestsTBB::PrintFuncTestName<TerekhovDGaussVertFuncTestsTBB>;
const auto kNameFnSTL = TerekhovDGaussVertFuncTestsSTL::PrintFuncTestName<TerekhovDGaussVertFuncTestsSTL>;

INSTANTIATE_TEST_SUITE_P(GaussFilterTestsSEQ, TerekhovDGaussVertFuncTestsSEQ, kSeqValues, kNameFnSEQ);
INSTANTIATE_TEST_SUITE_P(GaussFilterTestsOMP, TerekhovDGaussVertFuncTestsOMP, kOmpValues, kNameFnOMP);
INSTANTIATE_TEST_SUITE_P(GaussFilterTestsTBB, TerekhovDGaussVertFuncTestsTBB, kTbbValues, kNameFnTBB);
INSTANTIATE_TEST_SUITE_P(GaussFilterTestsSTL, TerekhovDGaussVertFuncTestsSTL, kStlValues, kNameFnSTL);

}  // namespace

}  // namespace terekhov_d_gauss_vert
