#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "krykov_e_sobel_op/common/include/common.hpp"
#include "krykov_e_sobel_op/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace krykov_e_sobel_op {

class KrykovERunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int test_id = std::get<0>(std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));

    const int size = 4;

    Image img;
    img.width = size;
    img.height = size;
    img.data.resize(static_cast<size_t>(size) * static_cast<size_t>(size));

    expected_output_.assign(static_cast<size_t>(size) * static_cast<size_t>(size), 0);

    switch (test_id) {
      case 0:
        SetUpConstantImage(img);
        break;
      case 1:
        SetUpVerticalEdge(img);
        break;
      case 2:
        SetUpHorizontalEdge(img);
        break;
      default:
        break;
    }

    input_data_ = img;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  static void SetUpConstantImage(Image &img) {
    for (auto &p : img.data) {
      p = {.r = 50, .g = 50, .b = 50};
    }
  }

  void SetUpVerticalEdge(Image &img) {
    const int size = img.width;
    for (int row = 0; row < size; ++row) {
      for (int col = 0; col < size; ++col) {
        uint8_t v = (col < 2) ? 0 : 255;
        img.data[(row * size) + col] = {.r = v, .g = v, .b = v};
      }
    }
    expected_output_[(1 * size) + 1] = 1020;
    expected_output_[(2 * size) + 1] = 1020;
    expected_output_[(1 * size) + 2] = 1020;
    expected_output_[(2 * size) + 2] = 1020;
  }

  void SetUpHorizontalEdge(Image &img) {
    const int size = img.width;
    for (int row = 0; row < size; ++row) {
      for (int col = 0; col < size; ++col) {
        uint8_t v = (row < 2) ? 0 : 255;
        img.data[(row * size) + col] = {.r = v, .g = v, .b = v};
      }
    }
    expected_output_[(1 * size) + 1] = 1020;
    expected_output_[(1 * size) + 2] = 1020;
    expected_output_[(2 * size) + 1] = 1020;
    expected_output_[(2 * size) + 2] = 1020;
  }

  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(KrykovERunFuncTestsThreads, SobelOp) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(0, "ConstantImage"), std::make_tuple(1, "VerticalEdge"),
                                            std::make_tuple(2, "HorizontalEdge")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<KrykovESobelOpSEQ, InType>(kTestParam, PPC_SETTINGS_krykov_e_sobel_op));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = KrykovERunFuncTestsThreads::PrintFuncTestName<KrykovERunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SobelTests, KrykovERunFuncTestsThreads, kGtestValues, kTestName);

}  // namespace
}  // namespace krykov_e_sobel_op
