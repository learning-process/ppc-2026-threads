#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
#include "chaschin_vladimir_linear_image_filtration_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace chaschin_v_linear_image_filtration_seq {

class ChaschinVRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  std::vector<float> static GenerateDeterministicImage(int width, int height) {
    std::vector<float> image(static_cast<std::vector<float>::size_type>(width) *
                             static_cast<std::vector<float>::size_type>(height));
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        image[(i * width) + j] = static_cast<float>((i + 1) * (j + 3) % 256);
      }
    }
    return image;
  }

  std::vector<float> ApplyGaussianKernel(const std::vector<float> &image, int width, int height) {
    std::vector<float> temp(
        static_cast<std::vector<float>::size_type>(width) * static_cast<std::vector<float>::size_type>(height), 0.0f);
    std::vector<float> output(
        static_cast<std::vector<float>::size_type>(width) * static_cast<std::vector<float>::size_type>(height), 0.0f);

    // Горизонтальный проход
    for (int yi = 0; yi < height; ++yi) {
      temp[(yi * width) + 0] = (image[(yi * width) + 0] * 2 + image[(yi * width) + 1]) / 3.F;
      for (int xy = 1; xy < width - 1; ++xy) {
        temp[(yi * width) + xy] =
            (image[(yi * width) + xy - 1] + 2.F * image[(yi * width) + xy] + image[(yi * width) + xy + 1]) / 4.F;
      }
      temp[(yi * width) + width - 1] = (image[(yi * width) + width - 2] + 2.F * image[(yi * width) + width - 1]) / 3.F;
    }

    // Вертикальный проход
    for (int xy = 0; xy < width; ++xy) {
      output[xy] = ((temp[xy] * 2) + temp[width + xy]) / 3.F;
      for (int yi = 1; yi < height - 1; ++yi) {
        output[(yi * width) + xy] =
            (temp[((yi - 1) * width) + xy] + 2.F * temp[(yi * width) + xy] + temp[((yi + 1) * width) + xy]) / 4.F;
      }
      output[((height - 1) * width) + xy] =
          (temp[((height - 2) * width) + xy] + 2.F * temp[((height - 1) * width) + xy]) / 3.F;
    }

    return output;
  }

  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);

    const int width = size;
    const int height = size;

    input_data_ = std::make_tuple(GenerateDeterministicImage(width, height), width, height);

    expected_output_ = ApplyGaussianKernel(std::get<0>(input_data_), width, height);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }

    constexpr float kEps = 1e-5F;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_output_[i]) > kEps) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(ChaschinVRunFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(4, "4"),   std::make_tuple(8, "8"),   std::make_tuple(16, "16"),
    std::make_tuple(32, "32"), std::make_tuple(64, "64"),
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<ChaschinVLinearFiltrationSEQ, InType>(kTestParam, PPC_SETTINGS_example_processes));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ChaschinVRunFuncTests::PrintFuncTestName<ChaschinVRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(LinearGaussianTests, ChaschinVRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chaschin_v_linear_image_filtration_seq
