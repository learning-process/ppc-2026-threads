#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "iskhakov_d_vertical_gauss_filter/common/include/common.hpp"
#include "iskhakov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"
#include "iskhakov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "iskhakov_d_vertical_gauss_filter/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace iskhakov_d_vertical_gauss_filter {

class IskhakovDVerticalGaussFilterFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<IskhakovDVerticalGaussFilterSEQ, InType>(
                                               kTestCases, PPC_SETTINGS_iskhakov_d_vertical_gauss_filter),
                                           ppc::util::AddFuncTask<IskhakovDVerticalGaussFilterOMP, InType>(
                                               kTestCases, PPC_SETTINGS_iskhakov_d_vertical_gauss_filter),
                                           ppc::util::AddFuncTask<IskhakovDVerticalGaussFilterTBB, InType>(
                                               kTestCases, PPC_SETTINGS_iskhakov_d_vertical_gauss_filter));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName =
    IskhakovDVerticalGaussFilterFuncTests::PrintFuncTestName<IskhakovDVerticalGaussFilterFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, IskhakovDVerticalGaussFilterFuncTests, kGtestValues, kFuncTestName);

}  // namespace

namespace {

struct InvalidInputParams {
  std::function<std::shared_ptr<BaseTask>(Matrix)> factory;
  Matrix input;
  std::string description;
};

void PrintTo(const InvalidInputParams &p, std::ostream *os) {
  *os << p.description;
}

class IskhakovDVerticalGaussFilterInvalidInputTests : public testing::TestWithParam<InvalidInputParams> {};

TEST_P(IskhakovDVerticalGaussFilterInvalidInputTests, Validation) {
  const auto &param = GetParam();
  auto task = param.factory(param.input);
  EXPECT_FALSE(task->Validation());
}

std::string PrintInvalidInputTestName(const testing::TestParamInfo<InvalidInputParams> &info) {
  return info.param.description;
}

Matrix zero_sizes{.width = 0, .height = 0, .data = {}};
Matrix zero_width{.width = 0, .height = 5, .data = std::vector<uint8_t>(5)};
Matrix zero_height{.width = 5, .height = 0, .data = std::vector<uint8_t>(5)};
Matrix size_mismatch{.width = 3, .height = 3, .data = {1, 2, 3}};
Matrix negative_width{.width = -1, .height = 5, .data = std::vector<uint8_t>(5)};
Matrix negative_height{.width = 5, .height = -1, .data = std::vector<uint8_t>(5)};

std::vector<InvalidInputParams> GenerateInvalidInputParams() {
  std::vector<InvalidInputParams> params;

  auto add = [&](const std::string &task_name, auto factory) {
    params.push_back({factory, zero_sizes, task_name + "_ZeroSizes"});
    params.push_back({factory, zero_width, task_name + "_ZeroWidthPositiveHeight"});
    params.push_back({factory, zero_height, task_name + "_ZeroHeightPositiveWidth"});
    params.push_back({factory, size_mismatch, task_name + "_DataSizeMismatch"});
    params.push_back({factory, negative_width, task_name + "_NegativeWidth"});
    params.push_back({factory, negative_height, task_name + "_NegativeHeight"});
  };

  add("SEQ", [](const Matrix &m) { return std::make_shared<IskhakovDVerticalGaussFilterSEQ>(m); });
  add("OMP", [](const Matrix &m) { return std::make_shared<IskhakovDVerticalGaussFilterOMP>(m); });
  add("TBB", [](const Matrix &m) { return std::make_shared<IskhakovDVerticalGaussFilterTBB>(m); });

  return params;
}

const auto kInvalidInputParams = GenerateInvalidInputParams();

INSTANTIATE_TEST_SUITE_P(InvalidInput, IskhakovDVerticalGaussFilterInvalidInputTests,
                         testing::ValuesIn(kInvalidInputParams), PrintInvalidInputTestName);

}  // namespace

}  // namespace iskhakov_d_vertical_gauss_filter
