#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "iskhakov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace iskhakov_d_vertical_gauss_filter {

class IskhakovDVerticalGaussFilterFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &input = std::get<0>(test_param);
    return std::to_string(input.width) + "x" + std::to_string(input.height);
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

TEST_P(IskhakovDVerticalGaussFilterFuncTests, RunTest) {
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

const std::array<TestType, 2> kTestCases = {std::make_tuple(kInput1x1, kExpected1x1),
                                            std::make_tuple(kInput2x2, kExpected2x2)};

using ParamType = std::tuple<std::function<std::shared_ptr<BaseTask>(InType)>, std::string, TestType>;

std::vector<ParamType> CreateTestParams() {
  std::vector<ParamType> params;
  for (const auto &test_case : kTestCases) {
    params.emplace_back([](const InType &in) -> std::shared_ptr<BaseTask> {
      return std::make_shared<IskhakovDVerticalGaussFilterSEQ>(in);
    }, "seq", test_case);
  }
  return params;
}

const auto kTestParams = CreateTestParams();
const auto kGtestValues = testing::ValuesIn(kTestParams);
const auto kFuncTestName =
    IskhakovDVerticalGaussFilterFuncTests::PrintFuncTestName<IskhakovDVerticalGaussFilterFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, IskhakovDVerticalGaussFilterFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace iskhakov_d_vertical_gauss_filter
