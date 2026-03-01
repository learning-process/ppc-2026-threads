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
  std::vector<float> GenerateDeterministicImage(int width, int height) {
    std::vector<float> image(static_cast<std::vector<float>::size_type>(width * height));
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        image[(i * width) + j] = static_cast<float>((i + 1) * (j + 3) % 256);
      }
    }
    return image;
  }

  std::vector<float> ApplyGaussianKernel(const std::vector<float> &image, int width, int height) {
    std::vector<float> temp(width * height, 0.0F);
    std::vector<float> output(width * height, 0.0F);

    // Горизонтальный проход
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        float left = (x > 0) ? image[y * width + (x - 1)] : image[y * width + x];
        float center = image[y * width + x];
        float right = (x < width - 1) ? image[y * width + (x + 1)] : image[y * width + x];
        temp[y * width + x] = (left + 2.F * center + right) / 4.F;
      }
    }

    // Вертикальный проход
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        float top = (y > 0) ? temp[(y - 1) * width + x] : temp[y * width + x];
        float center = temp[y * width + x];
        float bottom = (y < height - 1) ? temp[(y + 1) * width + x] : temp[y * width + x];
        output[y * width + x] = (top + 2.F * center + bottom) / 4.F;
      }
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
