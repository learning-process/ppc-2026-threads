#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "chaschin_vladimir_linear_image_filtration_seq/common/include/common.hpp"
#include "chaschin_vladimir_linear_image_filtration_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chaschin_v_linear_image_filtration_seq {

class ChaschinVRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kCount = 512;

  std::vector<float> static GenerateDeterministicImage(int width, int height) {
    std::vector<float> image(static_cast<std::vector<float>::size_type>(width * height));
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        image[(i * width) + j] = static_cast<float>((i + 1) * (j + 3) % 256);
      }
    }
    return image;
  }

  std::vector<float> static ApplyGaussianKernel(const std::vector<float> &image, int width, int height) {
    std::vector<float> temp(
        static_cast<std::vector<float>::size_type>(width) * static_cast<std::vector<float>::size_type>(height), 0.0F);
    std::vector<float> output
        static_cast<std::vector<float>::size_type>(width) * static_cast<std::vector<float>::size_type>(height), 0.0F);

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
    const int width = kCount;
    const int height = kCount;

    input_data_ = std::make_tuple(GenerateDeterministicImage(width, height), width, height);

    expected_output_ = ApplyGaussianKernel(std::get<0>(input_data_), width, height);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }

    constexpr float kEps = 1e-4F;
    for (size_t yi = 0; yi < output_data.size(); ++yi) {
      if (std::fabs(output_data[yi] - expected_output_[yi]) > kEps) {
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

TEST_P(ChaschinVRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ChaschinVLinearFiltrationSEQ>(PPC_SETTINGS_example_processes);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChaschinVRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ChaschinVRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chaschin_v_linear_image_filtration_seq
