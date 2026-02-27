#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>

#include "telnov_a_integral_rectangle/common/include/common.hpp"
#include "telnov_a_integral_rectangle/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace telnov_a_integral_rectangle {

class TelnovAIntegralRectangleFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &in = std::get<0>(test_param);
    return std::to_string(in.first) + "_" + std::to_string(in.second);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const int d = input_data_.second;
    const int n = input_data_.first;

    const double expected = static_cast<double>(d) / 2.0;
    const double tolerance = 1.0 / n;

    return std::abs(output_data - expected) < tolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(TelnovAIntegralRectangleFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    std::make_tuple(InType{1, 1}, "1D_N1"), std::make_tuple(InType{5, 1}, "1D_small"),
    std::make_tuple(InType{20, 1}, "1D"),   std::make_tuple(InType{10, 2}, "2D_small"),
    std::make_tuple(InType{20, 2}, "2D"),   std::make_tuple(InType{10, 3}, "3D_small"),
    std::make_tuple(InType{15, 3}, "3D"),   std::make_tuple(InType{8, 5}, "5D")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TelnovAIntegralRectangleSEQ, InType>(kTestParam, PPC_SETTINGS_telnov_a_integral_rectangle));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = TelnovAIntegralRectangleFuncTests::PrintFuncTestName<TelnovAIntegralRectangleFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, TelnovAIntegralRectangleFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace telnov_a_integral_rectangle
