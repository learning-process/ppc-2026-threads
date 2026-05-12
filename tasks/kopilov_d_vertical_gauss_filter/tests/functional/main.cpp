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

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "kopilov_d_vertical_gauss_filter/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class VerticalGaussFilterTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    const auto &input = std::get<0>(param);
    std::string data_suffix = input.data.empty() ? "Empty" : "Fill" + std::to_string(input.data[0]);
    return "Resolution_" + std::to_string(input.width) + "x" + std::to_string(input.height) + "_" + data_suffix;
  }

 protected:
  void SetUp() override {
    const auto &test_params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    in_matrix_ = std::get<0>(test_params);
    expected_matrix_ = std::get<1>(test_params);
  }

  bool CheckTestOutputData(OutType &actual_output) final {
    if (actual_output.width != expected_matrix_.width || actual_output.height != expected_matrix_.height) {
      return false;
    }
    return actual_output.data == expected_matrix_.data;
  }

  InType GetTestInputData() final {
    return in_matrix_;
  }

 private:
  InType in_matrix_;
  OutType expected_matrix_;
};

TEST_P(VerticalGaussFilterTest, CheckCorrectness) {
  ExecuteTest(GetParam());
}

namespace {

const std::array<TestType, 5> kValidationScenarios = {
    // 1x1
    std::make_tuple(Matrix{1, 1, {100}}, Matrix{1, 1, {100}}),
    // 2x2
    std::make_tuple(Matrix{2, 2, {1, 2, 3, 4}}, Matrix{2, 2, {1, 2, 2, 3}}),
    // 3x3
    std::make_tuple(Matrix{3, 3, std::vector<uint8_t>(9, 16)}, Matrix{3, 3, std::vector<uint8_t>(9, 16)}),
    std::make_tuple(Matrix{3, 3, std::vector<uint8_t>(9, 42)}, Matrix{3, 3, std::vector<uint8_t>(9, 42)}),
    // 4x4
    std::make_tuple(Matrix{4, 4, std::vector<uint8_t>(16, 100)}, Matrix{4, 4, std::vector<uint8_t>(16, 100)})};

const auto kTaskGenerators = std::tuple_cat(ppc::util::AddFuncTask<KopilovDVerticalGaussFilterSEQ, InType>(
                                                kValidationScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter),
                                            ppc::util::AddFuncTask<KopilovDVerticalGaussFilterOMP, InType>(
                                                kValidationScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter),
                                            ppc::util::AddFuncTask<KopilovDVerticalGaussFilterTBB, InType>(
                                                kValidationScenarios, PPC_SETTINGS_kopilov_d_vertical_gauss_filter));

INSTANTIATE_TEST_SUITE_P(FilterFunctionality, VerticalGaussFilterTest, ppc::util::ExpandToValues(kTaskGenerators),
                         VerticalGaussFilterTest::PrintFuncTestName<VerticalGaussFilterTest>);

}  // namespace

namespace {
void AssertValidationFails(const Matrix &m) {
  EXPECT_FALSE(std::make_shared<KopilovDVerticalGaussFilterSEQ>(m)->Validation());
  EXPECT_FALSE(std::make_shared<KopilovDVerticalGaussFilterOMP>(m)->Validation());
  EXPECT_FALSE(std::make_shared<KopilovDVerticalGaussFilterTBB>(m)->Validation());
}
}  // namespace

TEST(VerticalGaussFilterValidation, ZeroWidthFails) {
  AssertValidationFails(Matrix{0, 10, std::vector<uint8_t>(10, 0)});
}

TEST(VerticalGaussFilterValidation, ZeroHeightFails) {
  AssertValidationFails(Matrix{10, 0, std::vector<uint8_t>(10, 0)});
}

TEST(VerticalGaussFilterValidation, NegativeWidthFails) {
  AssertValidationFails(Matrix{-5, 5, std::vector<uint8_t>(25, 0)});
}

TEST(VerticalGaussFilterValidation, NegativeHeightFails) {
  AssertValidationFails(Matrix{5, -5, std::vector<uint8_t>(25, 0)});
}

TEST(VerticalGaussFilterValidation, DataSizeMismatchFails) {
  AssertValidationFails(Matrix{5, 5, std::vector<uint8_t>(15, 0)});
}

TEST(VerticalGaussFilterValidation, EmptyDataFails) {
  AssertValidationFails(Matrix{0, 0, {}});
}

}  // namespace kopilov_d_vertical_gauss_filter
