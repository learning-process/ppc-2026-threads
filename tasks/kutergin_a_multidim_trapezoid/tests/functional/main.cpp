#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"
#include "kutergin_a_multidim_trapezoid/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kutergin_a_multidim_trapezoid {
namespace {

auto f_const = [](const std::vector<double> &) { return 1.0; };

auto f_sum2 = [](const std::vector<double> &x) { return (x[0] + x[1]); };
auto f_sum3 = [](const std::vector<double> &x) { return (x[0] + x[1] + x[2]); };

auto f_square = [](const std::vector<double> &x) { return (x[0] * x[0]); };
auto f_square_sum = [](const std::vector<double> &x) { return ((x[0] * x[0]) + (x[1] * x[1])); };

auto f_mul2 = [](const std::vector<double> &x) { return (x[0] * x[1]); };

auto f_sin = [](const std::vector<double> &x) { return std::sin(x[0]); };
auto f_exp = [](const std::vector<double> &x) { return std::exp(x[0]); };

double Volume(const std::vector<std::pair<double, double>> &bounds) {
  double result = 1.0;
  for (const auto &[l, r] : bounds) {
    result *= (r - l);
  }
  return result;
}

double IntLinear(double a, double b) {
  return ((b * b - a * a) * 0.5);
}

double IntSquare(double a, double b) {
  return ((b * b * b - a * a * a) / 3.0);
}

double IntSin(double a, double b) {
  return (std::cos(a) - std::cos(b));
}

double IntExp(double a, double b) {
  return (std::exp(b) - std::exp(a));
}

class KuterginATrapezoidFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  using ParamType = typename ppc::util::BaseRunFuncTests<InType, OutType, TestType>::ParamType;

  static std::string GetName(const testing::TestParamInfo<ParamType> &info) {
    const auto &full_param = info.param;
    const auto &test_data = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(full_param);
    return std::get<0>(test_data);
  }

 protected:
  void SetUp() override {
    const auto &data = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    in_data_ = std::get<1>(data);
    ref_value_ = std::get<2>(data);
  }

  InType GetTestInputData() final {
    return in_data_;
  }

  bool CheckTestOutputData(OutType &output) final {
    constexpr double kTolerance = 1e-3;
    return (std::abs(output - ref_value_) <= kTolerance);
  }

 private:
  InType in_data_;
  double ref_value_ = 0.0;
};

TEST_P(KuterginATrapezoidFuncTest, Execution) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kFunctionalTests = {{
    {"Const_1D_unit", {f_const, {{0.0, 1.0}}, 120}, Volume({{0.0, 1.0}})},
    {"Const_2D_rect", {f_const, {{0.0, 2.0}, {1.0, 3.0}}, 100}, Volume({{0.0, 2.0}, {1.0, 3.0}})},
    {"Linear_1D", {[](const std::vector<double> &x) { return x[0]; }, {{-1.0, 2.0}}, 200}, IntLinear(-1.0, 2.0)},
    {"Linear_2D_sum", {f_sum2, {{0.0, 1.0}, {0.0, 1.0}}, 100}, (2.0 * IntLinear(0.0, 1.0) * 1.0)},
    {"Linear_3D_sum", {f_sum3, {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}, 60}, (3.0 * IntLinear(0.0, 1.0))},
    {"Square_1D", {f_square, {{0.0, 2.0}}, 200}, IntSquare(0.0, 2.0)},
    {"SquareSum_2D",
     {f_square_sum, {{0.0, 1.0}, {0.0, 2.0}}, 120},
     ((IntSquare(0.0, 1.0) * 2.0) + (IntSquare(0.0, 2.0) * 1.0))},
    {"Product_2D", {f_mul2, {{0.0, 2.0}, {0.0, 3.0}}, 150}, (IntLinear(0.0, 2.0) * IntLinear(0.0, 3.0))},
    {"Sin_0_pi", {f_sin, {{0.0, std::numbers::pi}}, 250}, IntSin(0.0, std::numbers::pi)},
    {"Exp_0_1", {f_exp, {{0.0, 1.0}}, 250}, IntExp(0.0, 1.0)},
    {"SinSin_2D",
     {[](const std::vector<double> &x) {
  return (std::sin(x[0]) * std::sin(x[1]));
}, {{0.0, std::numbers::pi}, {0.0, std::numbers::pi}}, 180},
     (IntSin(0.0, std::numbers::pi) * IntSin(0.0, std::numbers::pi))},
    {"Exp_2D",
     {[](const std::vector<double> &x) { return std::exp((x[0] + x[1])); }, {{0.0, 1.0}, {0.0, 1.0}}, 200},
     (std::pow(IntExp(0.0, 1.0), 2))},
}};

const auto kTaskPack = std::tuple_cat(ppc::util::AddFuncTask<KuterginAMultidimTrapezoidSEQ, InType>(
    kFunctionalTests, PPC_SETTINGS_kutergin_a_multidim_trapezoid));

const auto kValues = ppc::util::ExpandToValues(kTaskPack);

INSTANTIATE_TEST_SUITE_P(KuterginATrapezoidSuite, KuterginATrapezoidFuncTest, kValues,
                         KuterginATrapezoidFuncTest::GetName);

}  // namespace
}  // namespace kutergin_a_multidim_trapezoid
