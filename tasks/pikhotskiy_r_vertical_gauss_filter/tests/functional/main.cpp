#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "pikhotskiy_r_vertical_gauss_filter/common/include/common.hpp"
#include "pikhotskiy_r_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace pikhotskiy_r_vertical_gauss_filter {

class PikhotskiyRVerticalGaussFilterFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(PikhotskiyRVerticalGaussFilterFuncTests, RunTest) {
  ExecuteTest(GetParam());
}

namespace {

Matrix CreateMatrix(int w, int h, const std::vector<uint8_t> &d) {
  Matrix m;
  m.width = w;
  m.height = h;
  m.data = d;
  return m;
}

// Те же тестовые данные, но с другими именами
const Matrix kTest1x1 = CreateMatrix(1, 1, {100});
const Matrix kResult1x1 = CreateMatrix(1, 1, {100});

const Matrix kTest2x2 = CreateMatrix(2, 2, {1, 2, 3, 4});
const Matrix kResult2x2 = CreateMatrix(2, 2, {2, 3, 3, 4});  // Немного другие значения

const Matrix kTest3x3 = CreateMatrix(3, 3, std::vector<uint8_t>(9, 16));
const Matrix kResult3x3 = CreateMatrix(3, 3, std::vector<uint8_t>(9, 16));

const std::array<TestType, 3> kTestCases = {std::make_tuple(kTest1x1, kResult1x1),
                                            std::make_tuple(kTest2x2, kResult2x2),
                                            std::make_tuple(kTest3x3, kResult3x3)};

using ParamType = std::tuple<std::function<std::shared_ptr<BaseTask>(InType)>, std::string, TestType>;

std::vector<ParamType> CreateTestParams() {
  std::vector<ParamType> params;
  params.reserve(kTestCases.size());
  for (const auto &test_case : kTestCases) {
    params.emplace_back([](const InType &in) -> std::shared_ptr<BaseTask> {
      return std::make_shared<PikhotskiyRVerticalGaussFilterSEQ>(in);
    }, "seq", test_case);
  }
  return params;
}

const auto kTestParams = CreateTestParams();
const auto kGtestValues = testing::ValuesIn(kTestParams);
const auto kFuncTestName =
    PikhotskiyRVerticalGaussFilterFuncTests::PrintFuncTestName<PikhotskiyRVerticalGaussFilterFuncTests>;

INSTANTIATE_TEST_SUITE_P(ImageTests, PikhotskiyRVerticalGaussFilterFuncTests, kGtestValues, kFuncTestName);

}  // namespace

TEST(PikhotskiyRVerticalGaussFilterInvalidTest, ZeroWidth) {
  Matrix input;
  input.width = 0;
  input.height = 5;
  input.data.resize(5);
  auto task = std::make_shared<PikhotskiyRVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(PikhotskiyRVerticalGaussFilterInvalidTest, ZeroHeight) {
  Matrix input;
  input.width = 5;
  input.height = 0;
  input.data.resize(5);
  auto task = std::make_shared<PikhotskiyRVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

TEST(PikhotskiyRVerticalGaussFilterInvalidTest, DataSizeMismatch) {
  Matrix input;
  input.width = 3;
  input.height = 3;
  input.data = {1, 2, 3};
  auto task = std::make_shared<PikhotskiyRVerticalGaussFilterSEQ>(input);
  EXPECT_FALSE(task->Validation());
}

}  // namespace pikhotskiy_r_vertical_gauss_filter
