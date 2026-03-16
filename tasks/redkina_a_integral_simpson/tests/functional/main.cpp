#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"
#include "redkina_a_integral_simpson/omp/include/ops_omp.hpp"  // добавлено
#include "redkina_a_integral_simpson/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace redkina_a_integral_simpson {
namespace {

InputData MakeInput(double (*func)(const std::vector<double> &), std::vector<double> a, std::vector<double> b,
                    std::vector<int> n) {
  return InputData{.func = func, .a = std::move(a), .b = std::move(b), .n = std::move(n)};
}

TestType MakeTest(int id, double (*func)(const std::vector<double> &), std::vector<double> a, std::vector<double> b,
                  std::vector<int> n, double expected) {
  return std::make_tuple(id, func, std::move(a), std::move(b), std::move(n), expected);
}

class RedkinaAIntegralSimpsonFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "id_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = MakeInput(std::get<1>(params), std::get<2>(params), std::get<3>(params), std::get<4>(params));
    expected_ = std::get<5>(params);
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

const double kPi = std::numbers::pi;
const double kE = std::numbers::e;

const std::array<TestType, 20> kTestCases = {
    // 1D: константная функция
    MakeTest(1, [](const std::vector<double> &) { return 1.0; }, std::vector<double>{0.0}, std::vector<double>{1.0},
             std::vector<int>{2}, 1.0),
    // 1D: линейная функция
    MakeTest(2, [](const std::vector<double> &x) { return x[0]; }, std::vector<double>{0.0}, std::vector<double>{1.0},
             std::vector<int>{2}, 0.5),
    // 1D: квадратичная функция
    MakeTest(3, [](const std::vector<double> &x) { return x[0] * x[0]; }, std::vector<double>{0.0},
             std::vector<double>{1.0}, std::vector<int>{2}, 1.0 / 3.0),
    // 1D: кубическая функция
    MakeTest(4, [](const std::vector<double> &x) { return x[0] * x[0] * x[0]; }, std::vector<double>{0.0},
             std::vector<double>{1.0}, std::vector<int>{2}, 0.25),
    // 1D: x^4
    MakeTest(5, [](const std::vector<double> &x) { return x[0] * x[0] * x[0] * x[0]; }, std::vector<double>{0.0},
             std::vector<double>{1.0}, std::vector<int>{200}, 0.2),
    // 1D: sin(x)
    MakeTest(6, [](const std::vector<double> &x) { return std::sin(x[0]); }, std::vector<double>{0.0},
             std::vector<double>{kPi}, std::vector<int>{200}, 2.0),
    // 1D: exp(x)
    MakeTest(7, [](const std::vector<double> &x) { return std::exp(x[0]); }, std::vector<double>{0.0},
             std::vector<double>{1.0}, std::vector<int>{200}, kE - 1.0),
    // 2D: константная функция
    MakeTest(8, [](const std::vector<double> &) { return 1.0; }, std::vector<double>{0.0, 0.0},
             std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}, 1.0),
    // 2D: x*y
    MakeTest(9, [](const std::vector<double> &x) { return x[0] * x[1]; }, std::vector<double>{0.0, 0.0},
             std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}, 0.25),
    // 2D: x^2 + y
    MakeTest(10, [](const std::vector<double> &x) { return (x[0] * x[0]) + x[1]; }, std::vector<double>{0.0, 0.0},
             std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}, 5.0 / 6.0),
    // 2D: x*y^2
    MakeTest(11, [](const std::vector<double> &x) { return x[0] * (x[1] * x[1]); }, std::vector<double>{0.0, 0.0},
             std::vector<double>{1.0, 1.0}, std::vector<int>{2, 2}, 1.0 / 6.0),
    // 2D: exp(x+y)
    MakeTest(12, [](const std::vector<double> &x) { return std::exp(x[0] + x[1]); }, std::vector<double>{0.0, 0.0},
             std::vector<double>{1.0, 1.0}, std::vector<int>{200, 200}, (kE - 1.0) * (kE - 1.0)),
    // 2D: sin(x+y)
    MakeTest(13, [](const std::vector<double> &x) { return std::sin(x[0] + x[1]); }, std::vector<double>{0.0, 0.0},
             std::vector<double>{kPi, kPi}, std::vector<int>{200, 200}, 0.0),
    // 2D: sin(x)*cos(y)
    MakeTest(14, [](const std::vector<double> &x) { return std::sin(x[0]) * std::cos(x[1]); },
             std::vector<double>{0.0, 0.0}, std::vector<double>{kPi, kPi}, std::vector<int>{200, 200}, 0.0),
    // 2D: x*sin(y)
    MakeTest(15, [](const std::vector<double> &x) { return x[0] * std::sin(x[1]); }, std::vector<double>{0.0, 0.0},
             std::vector<double>{1.0, kPi}, std::vector<int>{200, 200}, 1.0),
    // 3D: константная функция
    MakeTest(16, [](const std::vector<double> &) { return 1.0; }, std::vector<double>{0.0, 0.0, 0.0},
             std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}, 1.0),
    // 3D: x*y*z
    MakeTest(17, [](const std::vector<double> &x) { return x[0] * x[1] * x[2]; }, std::vector<double>{0.0, 0.0, 0.0},
             std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}, 0.125),
    // 3D: x^2 + y^2 + z^2 на [0,1]^3
    MakeTest(18, [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); },
             std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}, 1.0),
    // 3D: x^2 + y^2 + z^2 на [-1,1]^3
    MakeTest(19, [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); },
             std::vector<double>{-1.0, -1.0, -1.0}, std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{2, 2, 2}, 8.0),
    // 3D: sin(x)*cos(y)*exp(z)
    MakeTest(20, [](const std::vector<double> &x) {
  return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]);
}, std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{kPi, kPi, 1.0}, std::vector<int>{40, 40, 40}, 0.0)};

// Для последовательной версии
const auto kTestTasksListSeq =
    ppc::util::AddFuncTask<RedkinaAIntegralSimpsonSEQ, InType>(kTestCases, PPC_SETTINGS_redkina_a_integral_simpson);
const auto kGtestValuesSeq = ppc::util::ExpandToValues(kTestTasksListSeq);

// Для OMP версии
const auto kTestTasksListOmp =
    ppc::util::AddFuncTask<RedkinaAIntegralSimpsonOMP, InType>(kTestCases, PPC_SETTINGS_redkina_a_integral_simpson);
const auto kGtestValuesOmp = ppc::util::ExpandToValues(kTestTasksListOmp);

const auto kTestName = RedkinaAIntegralSimpsonFuncTests::PrintFuncTestName<RedkinaAIntegralSimpsonFuncTests>;

INSTANTIATE_TEST_SUITE_P(IntegralSimpsonTestsSeq, RedkinaAIntegralSimpsonFuncTests, kGtestValuesSeq, kTestName);
INSTANTIATE_TEST_SUITE_P(IntegralSimpsonTestsOmp, RedkinaAIntegralSimpsonFuncTests, kGtestValuesOmp, kTestName);

TEST_P(RedkinaAIntegralSimpsonFuncTests, CheckIntegral) {
  ExecuteTest(GetParam());
}

}  // namespace
}  // namespace redkina_a_integral_simpson
