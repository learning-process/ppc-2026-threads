#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "redkina_a_integral_simpson_seq/common/include/common.hpp"
#include "redkina_a_integral_simpson_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace redkina_a_integral_simpson_seq {

class RedkinaAIntegralSimpsonFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "id_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<1>(params);
    expected_ = std::get<2>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-6;
    return std::fabs(output_data - expected_) < eps;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  double expected_ = 0.0;
};

namespace {

const double kPi = std::numbers::pi;
const double kE = std::numbers::e;

InputData MakeInput(std::function<double(const std::vector<double> &)> func, std::vector<double> a,
                    std::vector<double> b, std::vector<int> n) {
  return InputData{.func = std::move(func), .a = std::move(a), .b = std::move(b), .n = std::move(n)};
}

const std::array<TestType, 20> kTestCases = {
    {// 1D: константа
     std::make_tuple(1,
                     MakeInput([](const std::vector<double> &) { return 1.0; }, std::vector<double>{0.0},
                               std::vector<double>{1.0}, std::vector<int>{2}),
                     1.0),

     // 1D: линейная
     std::make_tuple(2,
                     MakeInput([](const std::vector<double> &x) { return x[0]; }, std::vector<double>{0.0},
                               std::vector<double>{1.0}, std::vector<int>{2}),
                     0.5),

     // 1D: квадратичная
     std::make_tuple(3,
                     MakeInput([](const std::vector<double> &x) { return x[0] * x[0]; }, std::vector<double>{0.0},
                               std::vector<double>{1.0}, std::vector<int>{2}),
                     1.0 / 3.0),

     // 1D: кубическая
     std::make_tuple(4,
                     MakeInput([](const std::vector<double> &x) { return x[0] * x[0] * x[0]; },
                               std::vector<double>{0.0}, std::vector<double>{1.0}, std::vector<int>{2}),
                     0.25),

     // 1D: x^4
     std::make_tuple(5,
                     MakeInput([](const std::vector<double> &x) { return x[0] * x[0] * x[0] * x[0]; },
                               std::vector<double>{0.0}, std::vector<double>{1.0}, std::vector<int>{200}),
                     0.2),

     // 1D: sin(x)
     std::make_tuple(6,
                     MakeInput([](const std::vector<double> &x) { return std::sin(x[0]); }, std::vector<double>{0.0},
                               std::vector<double>{kPi}, std::vector<int>{200}),
                     2.0),

     // 1D: exp(x)
     std::make_tuple(7,
                     MakeInput([](const std::vector<double> &x) { return std::exp(x[0]); }, std::vector<double>{0.0},
                               std::vector<double>{1.0}, std::vector<int>{200}),
                     kE - 1.0),

     // 2D: константа
     std::make_tuple(8,
                     MakeInput([](const std::vector<double> &) { return 1.0; }, std::vector<double>{0.0, 0.0},
                               std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}),
                     1.0),

     // 2D: x*y
     std::make_tuple(9,
                     MakeInput([](const std::vector<double> &x) { return x[0] * x[1]; }, std::vector<double>{0.0, 0.0},
                               std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}),
                     0.25),

     // 2D: x^2 + y
     std::make_tuple(10,
                     MakeInput([](const std::vector<double> &x) { return (x[0] * x[0]) + x[1]; },
                               std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}),
                     5.0 / 6.0),

     // 2D: x*y^2
     std::make_tuple(11,
                     MakeInput([](const std::vector<double> &x) { return x[0] * (x[1] * x[1]); },
                               std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}),
                     1.0 / 6.0),

     // 2D: exp(x+y)
     std::make_tuple(
         12,
         MakeInput([](const std::vector<double> &x) { return std::exp(x[0] + x[1]); }, std::vector<double>{0.0, 0.0},
                   std::vector<double>{1.0, 1.0}, std::vector<int>{200, 200}),
         (kE - 1.0) * (kE - 1.0)),

     // 2D: sin(x+y)
     std::make_tuple(
         13,
         MakeInput([](const std::vector<double> &x) { return std::sin(x[0] + x[1]); }, std::vector<double>{0.0, 0.0},
                   std::vector<double>{kPi, kPi}, std::vector<int>{200, 200}),
         0.0),

     // 2D: sin(x)*cos(y)
     std::make_tuple(
         14,
         MakeInput([](const std::vector<double> &x) { return std::sin(x[0]) * std::cos(x[1]); },
                   std::vector<double>{0.0, 0.0}, std::vector<double>{kPi, kPi}, std::vector<int>{200, 200}),
         0.0),

     // 2D: x*sin(y)
     std::make_tuple(
         15,
         MakeInput([](const std::vector<double> &x) { return x[0] * std::sin(x[1]); }, std::vector<double>{0.0, 0.0},
                   std::vector<double>{1.0, kPi}, std::vector<int>{200, 200}),
         1.0),

     // 3D: константа
     std::make_tuple(16,
                     MakeInput([](const std::vector<double> &) { return 1.0; }, std::vector<double>{0.0, 0.0, 0.0},
                               std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}),
                     1.0),

     // 3D: x*y*z
     std::make_tuple(
         17,
         MakeInput([](const std::vector<double> &x) { return x[0] * x[1] * x[2]; }, std::vector<double>{0.0, 0.0, 0.0},
                   std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}),
         0.125),

     // 3D: x^2 + y^2 + z^2 на [0,1]^3
     std::make_tuple(
         18,
         MakeInput([](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); },
                   std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}),
         1.0),

     // 3D: x^2 + y^2 + z^2 на [-1,1]^3
     std::make_tuple(
         19,
         MakeInput([](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); },
                   std::vector<double>{-1.0, -1.0, -1.0}, std::vector<double>{1.0, 1.0, 1.0},
                   std::vector<int>{2, 2, 2}),
         8.0),

     // 3D: sin(x)*cos(y)*exp(z) (уменьшенное разбиение для скорости)
     std::make_tuple(
         20,
         MakeInput([](const std::vector<double> &x) { return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]); },
                   std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{kPi, kPi, 1.0},
                   std::vector<int>{40, 40, 40}),
         0.0)}};

const auto kTestTasksList =
    ppc::util::AddFuncTask<RedkinaAIntegralSimpsonSEQ, InType>(kTestCases, PPC_SETTINGS_redkina_a_integral_simpson_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = RedkinaAIntegralSimpsonFuncTests::PrintFuncTestName<RedkinaAIntegralSimpsonFuncTests>;

INSTANTIATE_TEST_SUITE_P(IntegralSimpsonTests, RedkinaAIntegralSimpsonFuncTests, kGtestValues, kTestName);

TEST_P(RedkinaAIntegralSimpsonFuncTests, CheckIntegral) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace redkina_a_integral_simpson_seq
