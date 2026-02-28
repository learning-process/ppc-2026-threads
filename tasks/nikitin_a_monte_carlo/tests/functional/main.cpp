#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "nikitin_a_monte_carlo/common/include/common.hpp"
#include "nikitin_a_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nikitin_a_monte_carlo {

class NikitinAMonteCarloFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    auto params = GetParam();
    TestType test_params = std::get<2>(params);

    test_case_id_ = std::get<0>(test_params);
    test_description_ = std::get<1>(test_params);
    input_data_ = std::get<2>(test_params);

    // Вычисляем точные значения интегралов для разных тестов
    switch (test_case_id_) {
      // 1D тесты (1-5)
      case 1:
        expected_output_ = 1.0;
        tolerance_ = 0.02;
        break;
      case 2:
        expected_output_ = 10.0;
        tolerance_ = 0.02;
        break;
      case 3:
        expected_output_ = 0.5;
        tolerance_ = 0.02;
        break;
      case 4:
        expected_output_ = 2.0;
        tolerance_ = 0.02;
        break;
      case 5:
        expected_output_ = std::exp(1.0) - 1.0;
        tolerance_ = 0.03;
        break;

      // 2D тесты (6-12)
      case 6:
        expected_output_ = 1.0;
        tolerance_ = 0.02;
        break;
      case 7:
        expected_output_ = 6.0;
        tolerance_ = 0.02;
        break;
      case 8:
        expected_output_ = 0.5;
        tolerance_ = 0.03;
        break;
      case 9:
        expected_output_ = 0.25;
        tolerance_ = 0.03;
        break;
      case 10:
        expected_output_ = 9.0;
        tolerance_ = 0.03;
        break;
      case 11:
        expected_output_ = 2.0 / 3.0;
        tolerance_ = 0.03;
        break;
      case 12:
        expected_output_ = 8.0 / 3.0;
        tolerance_ = 0.03;
        break;

      // 3D тесты (13-15)
      case 13:
        expected_output_ = 1.0;
        tolerance_ = 0.02;
        break;
      case 14:
        expected_output_ = 24.0;
        tolerance_ = 0.02;
        break;
      case 15:
        expected_output_ = 0.5;
        tolerance_ = 0.03;
        break;

      // Тесты с разным количеством точек (16-18)
      case 16:
        expected_output_ = 1.0;
        tolerance_ = 0.2;
        break;
      case 17:
        expected_output_ = 1.0;
        tolerance_ = 0.05;
        break;
      case 18:
        expected_output_ = 1.0;
        tolerance_ = 0.01;
        break;

      // Комбинированные тесты (19-22)
      case 19:
        expected_output_ = 10.5;
        tolerance_ = 0.03;
        break;
      case 20:
        expected_output_ = 3.0;
        tolerance_ = 0.02;
        break;
      case 21:
        expected_output_ = -4.0;
        tolerance_ = 0.03;
        break;
      case 22:
        expected_output_ = 42.0;
        tolerance_ = 0.04;
        break;

      default:
        throw std::runtime_error("Unknown test case ID: " + std::to_string(test_case_id_));
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (std::abs(expected_output_) > 1e-10) {
      double relative_error = std::abs(output_data - expected_output_) / std::abs(expected_output_);
      return relative_error <= tolerance_;
    } else {
      return std::abs(output_data - expected_output_) <= tolerance_;
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  int test_case_id_ = 0;
  std::string test_description_;
  InType input_data_;
  double expected_output_ = 0.0;
  double tolerance_ = 0.01;
};

namespace {

TEST_P(NikitinAMonteCarloFuncTests, MonteCarloIntegrationTest) {
  ExecuteTest(GetParam());
}

// Только корректные тесты (без тестов валидации)
const std::array<TestType, 22> kTestParam = {
    {// 1D тесты (1-5)
     TestType{1, "1d_constant_0_1",
              std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, 10000, FunctionType::kConstant)},
     TestType{2, "1d_constant_0_10",
              std::make_tuple(std::vector<double>{0.0}, std::vector<double>{10.0}, 10000, FunctionType::kConstant)},
     TestType{3, "1d_linear_0_1",
              std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, 10000, FunctionType::kLinear)},
     TestType{4, "1d_linear_0_2",
              std::make_tuple(std::vector<double>{0.0}, std::vector<double>{2.0}, 10000, FunctionType::kLinear)},
     TestType{5, "1d_exponential_0_1",
              std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, 20000, FunctionType::kExponential)},

     // 2D тесты (6-12)
     TestType{
         6, "2d_constant_unit",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 10000, FunctionType::kConstant)},
     TestType{
         7, "2d_constant_rectangle",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{2.0, 3.0}, 10000, FunctionType::kConstant)},
     TestType{
         8, "2d_linear_unit",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 10000, FunctionType::kLinear)},
     TestType{
         9, "2d_product_unit",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 15000, FunctionType::kProduct)},
     TestType{
         10, "2d_product_rectangle",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{2.0, 3.0}, 15000, FunctionType::kProduct)},
     TestType{11, "2d_quadratic_unit",
              std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 15000,
                              FunctionType::kQuadratic)},
     TestType{12, "2d_quadratic_symmetric",
              std::make_tuple(std::vector<double>{-1.0, -1.0}, std::vector<double>{1.0, 1.0}, 15000,
                              FunctionType::kQuadratic)},

     // 3D тесты (13-15)
     TestType{13, "3d_constant_unit",
              std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, 20000,
                              FunctionType::kConstant)},
     TestType{14, "3d_constant_rectangular",
              std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{2.0, 3.0, 4.0}, 20000,
                              FunctionType::kConstant)},
     TestType{15, "3d_linear_unit",
              std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, 20000,
                              FunctionType::kLinear)},

     // Тесты с разным количеством точек (16-18)
     TestType{
         16, "pts_few_100",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 100, FunctionType::kConstant)},
     TestType{
         17, "pts_medium_1000",
         std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 1000, FunctionType::kConstant)},
     TestType{18, "pts_many_100000",
              std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 100000,
                              FunctionType::kConstant)},

     // Комбинированные тесты (19-22)
     TestType{19, "comb_shifted_linear",
              std::make_tuple(std::vector<double>{2.0}, std::vector<double>{5.0}, 10000, FunctionType::kLinear)},
     TestType{20, "comb_negative_constant",
              std::make_tuple(std::vector<double>{-5.0}, std::vector<double>{-2.0}, 10000, FunctionType::kConstant)},
     TestType{21, "comb_negative_linear",
              std::make_tuple(std::vector<double>{-3.0}, std::vector<double>{-1.0}, 10000, FunctionType::kLinear)},
     TestType{22, "comb_asymmetric_product",
              std::make_tuple(std::vector<double>{1.0, 2.0}, std::vector<double>{3.0, 5.0}, 20000,
                              FunctionType::kProduct)}}};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<NikitinAMonteCarloSEQ, InType>(kTestParam, PPC_SETTINGS_nikitin_a_monte_carlo));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = NikitinAMonteCarloFuncTests::PrintFuncTestName<NikitinAMonteCarloFuncTests>;

INSTANTIATE_TEST_SUITE_P(MonteCarloTests, NikitinAMonteCarloFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nikitin_a_monte_carlo
