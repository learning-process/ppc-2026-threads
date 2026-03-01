#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class KopilovDVerticalGaussFilterFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &input = std::get<0>(test_param);
    std::string extra = input.data.empty() ? "empty" : std::to_string(input.data[0]);
    return std::to_string(input.width) + "x" + std::to_string(input.height) + "_" + extra;
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_data_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.width == expected_data_.width && output_data.height == expected_data_.height &&
           output_data.data == expected_data_.data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

TEST_P(KopilovDVerticalGaussFilterFuncTests, RunTest) {
  ExecuteTest(GetParam());
}

namespace {

Matrix MakeMatrix(int width, int height, const std::vector<uint8_t> &data) {
  Matrix m;
  m.width = width;
  m.height = height;
  m.data = data;
  return m;
}

const Matrix kInput1x1 = MakeMatrix(1, 1, {100});
const Matrix kExpected1x1 = MakeMatrix(1, 1, {100});

const Matrix kInput2x2 = MakeMatrix(2, 2, {1, 2, 3, 4});
const Matrix kExpected2x2 = MakeMatrix(2, 2, {1, 2, 2, 3});

const Matrix kInput3x316 = MakeMatrix(3, 3, std::vector<uint8_t>(9, 16));
const Matrix kExpected3x316 = MakeMatrix(3, 3, std::vector<uint8_t>(9, 16));

const Matrix kInput3x342 = MakeMatrix(3, 3, std::vector<uint8_t>(9, 42));
const Matrix kExpected3x342 = MakeMatrix(3, 3, std::vector<uint8_t>(9, 42));

const Matrix kInput4x4100 = MakeMatrix(4, 4, std::vector<uint8_t>(16, 100));
const Matrix kExpected4x4100 = MakeMatrix(4, 4, std::vector<uint8_t>(16, 100));

const std::array<TestType, 5> kTestCases = {
    std::make_tuple(kInput1x1, kExpected1x1), std::make_tuple(kInput2x2, kExpected2x2),
    std::make_tuple(kInput3x316, kExpected3x316), std::make_tuple(kInput3x342, kExpected3x342),
    std::make_tuple(kInput4x4100, kExpected4x4100)};

using ParamType = std::tuple<std::function<std::shared_ptr<BaseTask>(InType)>, std::string, TestType>;

std::vector<ParamType> CreateTestParams() {
  std::vector<ParamType> params;
  params.reserve(kTestCases.size());
  for (const auto &test_case : kTestCases) {
    params.emplace_back([](const InType &in) -> std::shared_ptr<BaseTask> {
      return std::make_shared<KopilovDVerticalGaussFilterSEQ>(in);
    }, "seq", test_case);
  }
  return params;
}

const auto kTestParams = CreateTestParams();
const auto kGtestValues = testing::ValuesIn(kTestParams);
const auto kFuncTestName =
    KopilovDVerticalGaussFilterFuncTests::PrintFuncTestName<KopilovDVerticalGaussFilterFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, KopilovDVerticalGaussFilterFuncTests, kGtestValues, kFuncTestName);

}  // namespace

TEST(KopilovDVerticalGaussFilterInvalidInputTest, ZeroSizes) {
  Matrix input;
  input.width = 0;
  input.height = 0;
  input.data = {};
  auto task = std::make_shared<KopilovDVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(KopilovDVerticalGaussFilterInvalidInputTest, ZeroWidthPositiveHeight) {
  Matrix input;
  input.width = 0;
  input.height = 5;
  input.data.resize(5);
  auto task = std::make_shared<KopilovDVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(KopilovDVerticalGaussFilterInvalidInputTest, ZeroHeightPositiveWidth) {
  Matrix input;
  input.width = 5;
  input.height = 0;
  input.data.resize(5);
  auto task = std::make_shared<KopilovDVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(KopilovDVerticalGaussFilterInvalidInputTest, DataSizeMismatch) {
  Matrix input;
  input.width = 3;
  input.height = 3;
  input.data = {1, 2, 3};
  auto task = std::make_shared<KopilovDVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(KopilovDVerticalGaussFilterInvalidInputTest, NegativeWidth) {
  Matrix input;
  input.width = -1;
  input.height = 5;
  input.data.resize(5);
  auto task = std::make_shared<KopilovDVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(KopilovDVerticalGaussFilterInvalidInputTest, NegativeHeight) {
  Matrix input;
  input.width = 5;
  input.height = -1;
  input.data.resize(5);
  auto task = std::make_shared<KopilovDVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

}  // namespace kopilov_d_vertical_gauss_filter
