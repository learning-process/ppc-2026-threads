#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "balchunayte_z_sobel/common/include/common.hpp"
#include "balchunayte_z_sobel/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace balchunayte_z_sobel {

class BalchunayteZRunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const int test_id = std::get<0>(std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    const int image_size = 4;

    Image input_image;
    input_image.width = image_size;
    input_image.height = image_size;
    input_image.data.resize(static_cast<size_t>(image_size) * static_cast<size_t>(image_size));

    expected_output_.assign(static_cast<size_t>(image_size) * static_cast<size_t>(image_size), 0);

    if (test_id == 0) {
      SetUpConstantImage(input_image);
    } else if (test_id == 1) {
      SetUpVerticalEdge(input_image);
    } else if (test_id == 2) {
      SetUpHorizontalEdge(input_image);
    }

    input_data_ = input_image;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  static void SetUpConstantImage(Image &input_image) {
    for (auto &pixel_value : input_image.data) {
      pixel_value = Pixel{50, 50, 50};
    }
  }

  void SetUpVerticalEdge(Image &input_image) {
    const int image_size = input_image.width;

    for (int row_index = 0; row_index < image_size; ++row_index) {
      for (int col_index = 0; col_index < image_size; ++col_index) {
        const uint8_t intensity_value = (col_index < 2) ? 0 : 255;
        input_image.data[static_cast<size_t>((row_index * image_size) + col_index)] =
            Pixel{intensity_value, intensity_value, intensity_value};
      }
    }

    expected_output_[static_cast<size_t>((1 * image_size) + 1)] = 1020;
    expected_output_[static_cast<size_t>((2 * image_size) + 1)] = 1020;
    expected_output_[static_cast<size_t>((1 * image_size) + 2)] = 1020;
    expected_output_[static_cast<size_t>((2 * image_size) + 2)] = 1020;
  }

  void SetUpHorizontalEdge(Image &input_image) {
    const int image_size = input_image.width;

    for (int row_index = 0; row_index < image_size; ++row_index) {
      for (int col_index = 0; col_index < image_size; ++col_index) {
        const uint8_t intensity_value = (row_index < 2) ? 0 : 255;
        input_image.data[static_cast<size_t>((row_index * image_size) + col_index)] =
            Pixel{intensity_value, intensity_value, intensity_value};
      }
    }

    expected_output_[static_cast<size_t>((1 * image_size) + 1)] = 1020;
    expected_output_[static_cast<size_t>((1 * image_size) + 2)] = 1020;
    expected_output_[static_cast<size_t>((2 * image_size) + 1)] = 1020;
    expected_output_[static_cast<size_t>((2 * image_size) + 2)] = 1020;
  }

  InType input_data_{};
  OutType expected_output_{};
};

namespace {

TEST_P(BalchunayteZRunFuncTestsSEQ, SobelOp) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(0, "ConstantImage"), std::make_tuple(1, "VerticalEdge"),
                                            std::make_tuple(2, "HorizontalEdge")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BalchunayteZSobelOpSEQ, InType>(kTestParam, PPC_SETTINGS_balchunayte_z_sobel));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = BalchunayteZRunFuncTestsSEQ::PrintFuncTestName<BalchunayteZRunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(SobelTests, BalchunayteZRunFuncTestsSEQ, kGtestValues, kTestName);

}  // namespace

}  // namespace balchunayte_z_sobel
