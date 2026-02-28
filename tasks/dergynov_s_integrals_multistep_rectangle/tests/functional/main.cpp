#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"
#include "dergynov_s_integrals_multistep_rectangle/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {
auto const_func = [](const std::vector<double> &) { return 1.0; };

auto linear_1d = [](const std::vector<double> &x) { return x[0]; };
auto linear_2d = [](const std::vector<double> &x) { return x[0] + x[1]; };
auto linear_3d = [](const std::vector<double> &x) { return x[0] + x[1] + x[2]; };

auto product_2d = [](const std::vector<double> &x) { return x[0] * x[1]; };
auto product_3d = [](const std::vector<double> &x) { return x[0] * x[1] * x[2]; };

auto quadratic_1d = [](const std::vector<double> &x) { return x[0] * x[0]; };
auto quadratic_2d = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]); };

auto sin_1d = [](const std::vector<double> &x) { return std::sin(x[0]); };
auto cos_1d = [](const std::vector<double> &x) { return std::cos(x[0]); };
auto sin_sin_2d = [](const std::vector<double> &x) { return std::sin(x[0]) * std::sin(x[1]); };

auto exp_1d = [](const std::vector<double> &x) { return std::exp(x[0]); };
auto exp_2d = [](const std::vector<double> &x) { return std::exp(x[0] + x[1]); };

double ExactConst(const std::vector<std::pair<double, double>> &borders) {
  double vol = 1.0;
  for (const auto &[a, b] : borders) {
    vol *= (b - a);
  }
  return vol;
}

double ExactLinear1D(double a, double b) {
  return (b * b - a * a) / 2.0;
}

double ExactLinear2D(double a, double b) {
  double i1 = (b * b - a * a) / 2.0;
  double len = b - a;
  return 2.0 * i1 * len;
}

double ExactLinear3D(double a, double b) {
  double i1 = (b * b - a * a) / 2.0;
  double len = b - a;
  return 3.0 * i1 * len * len;
}

double ExactProduct2D(double a1, double b1, double a2, double b2) {
  double i1 = (b1 * b1 - a1 * a1) / 2.0;
  double i2 = (b2 * b2 - a2 * a2) / 2.0;
  return i1 * i2;
}

double ExactProduct3D(double a1, double b1, double a2, double b2, double a3, double b3) {
  double i1 = (b1 * b1 - a1 * a1) / 2.0;
  double i2 = (b2 * b2 - a2 * a2) / 2.0;
  double i3 = (b3 * b3 - a3 * a3) / 2.0;
  return i1 * i2 * i3;
}

double ExactQuadratic1D(double a, double b) {
  return (b * b * b - a * a * a) / 3.0;
}

double ExactQuadratic2D(double a1, double b1, double a2, double b2) {
  double int_x2 = (b1 * b1 * b1 - a1 * a1 * a1) / 3.0;
  double int_y2 = (b2 * b2 * b2 - a2 * a2 * a2) / 3.0;
  double len_x = b1 - a1;
  double len_y = b2 - a2;
  return (int_x2 * len_y) + (int_y2 * len_x);
}

double ExactSin1D(double a, double b) {
  return std::cos(a) - std::cos(b);
}

double ExactSinSin2D(double a1, double b1, double a2, double b2) {
  return ExactSin1D(a1, b1) * ExactSin1D(a2, b2);
}

double ExactExp1D(double a, double b) {
  return std::exp(b) - std::exp(a);
}

double ExactExp2D(double a, double b) {
  double i1 = std::exp(b) - std::exp(a);
  return i1 * i1;
}

class DergynovSIntegralsFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<1>(params);
    expected_ = std::get<2>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-3;
    return std::fabs(output_data - expected_) <= eps;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_ = 0.0;
};

TEST_P(DergynovSIntegralsFuncTest, Run) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 19> kTests = {{
    TestType{"Const_1D_0_1_100", InType{const_func, std::vector<std::pair<double, double>>{{0.0, 1.0}}, 100},
             ExactConst({{0.0, 1.0}})},

    TestType{"Const_1D_m2_3_200", InType{const_func, std::vector<std::pair<double, double>>{{-2.0, 3.0}}, 200},
             ExactConst({{-2.0, 3.0}})},

    TestType{"Const_2D_0_1_x_0_2_80",
             InType{const_func, std::vector<std::pair<double, double>>{{0.0, 1.0}, {0.0, 2.0}}, 80},
             ExactConst({{0.0, 1.0}, {0.0, 2.0}})},

    TestType{"Const_3D_0_1_x_m1_1_x_2_5_40",
             InType{const_func, std::vector<std::pair<double, double>>{{0.0, 1.0}, {-1.0, 1.0}, {2.0, 5.0}}, 40},
             ExactConst({{0.0, 1.0}, {-1.0, 1.0}, {2.0, 5.0}})},

    TestType{"Linear_1D_0_2_200", InType{linear_1d, std::vector<std::pair<double, double>>{{0.0, 2.0}}, 200},
             ExactLinear1D(0.0, 2.0)},

    TestType{"Linear_1D_m1_1_200", InType{linear_1d, std::vector<std::pair<double, double>>{{-1.0, 1.0}}, 200},
             ExactLinear1D(-1.0, 1.0)},

    TestType{"Linear_2D_x_plus_y_0_1_80",
             InType{linear_2d, std::vector<std::pair<double, double>>{{0.0, 1.0}, {0.0, 1.0}}, 80},
             ExactLinear2D(0.0, 1.0)},

    TestType{"Linear_2D_x_plus_y_m1_1_80",
             InType{linear_2d, std::vector<std::pair<double, double>>{{-1.0, 1.0}, {-1.0, 1.0}}, 80},
             ExactLinear2D(-1.0, 1.0)},

    TestType{"Linear_3D_sum_0_1_40",
             InType{linear_3d, std::vector<std::pair<double, double>>{{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}, 40},
             ExactLinear3D(0.0, 1.0)},

    TestType{"Product_2D_xy_0_2_x_1_3_100",
             InType{product_2d, std::vector<std::pair<double, double>>{{0.0, 2.0}, {1.0, 3.0}}, 100},
             ExactProduct2D(0.0, 2.0, 1.0, 3.0)},

    TestType{"Product_2D_xy_m1_1_x_0_2_100",
             InType{product_2d, std::vector<std::pair<double, double>>{{-1.0, 1.0}, {0.0, 2.0}}, 100}, 0.0},

    TestType{"Product_3D_xyz_0_1_x_0_2_x_0_3_100",
             InType{product_3d, std::vector<std::pair<double, double>>{{0.0, 1.0}, {0.0, 2.0}, {0.0, 3.0}}, 100},
             ExactProduct3D(0.0, 1.0, 0.0, 2.0, 0.0, 3.0)},

    TestType{"Quadratic_1D_x2_m1_1_200", InType{quadratic_1d, std::vector<std::pair<double, double>>{{-1.0, 1.0}}, 200},
             ExactQuadratic1D(-1.0, 1.0)},

    TestType{"Quadratic_2D_x2_plus_y2_0_1_x_0_2_100",
             InType{quadratic_2d, std::vector<std::pair<double, double>>{{0.0, 1.0}, {0.0, 2.0}}, 100},
             ExactQuadratic2D(0.0, 1.0, 0.0, 2.0)},

    TestType{"Sin_1D_0_pi_200", InType{sin_1d, std::vector<std::pair<double, double>>{{0.0, std::numbers::pi}}, 200},
             ExactSin1D(0.0, std::numbers::pi)},

    TestType{"Cos_1D_0_pi2_200",
             InType{cos_1d, std::vector<std::pair<double, double>>{{0.0, std::numbers::pi / 2.0}}, 200},
             ExactSin1D(0.0, std::numbers::pi / 2.0)},

    TestType{
        "SinSin_2D_0_pi_x_0_pi2_200",
        InType{sin_sin_2d,
               std::vector<std::pair<double, double>>{{0.0, std::numbers::pi}, {0.0, std::numbers::pi / 2.0}}, 200},
        ExactSinSin2D(0.0, std::numbers::pi, 0.0, std::numbers::pi / 2.0)},

    TestType{"Exp_1D_0_1_200", InType{exp_1d, std::vector<std::pair<double, double>>{{0.0, 1.0}}, 200},
             ExactExp1D(0.0, 1.0)},

    TestType{"Exp_2D_0_1_200", InType{exp_2d, std::vector<std::pair<double, double>>{{0.0, 1.0}, {0.0, 1.0}}, 200},
             ExactExp2D(0.0, 1.0)},
}};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<DergynovSIntegralsMultistepRectangleSEQ, InType>(
    kTests, PPC_SETTINGS_dergynov_s_integrals_multistep_rectangle));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = DergynovSIntegralsFuncTest::PrintFuncTestName<DergynovSIntegralsFuncTest>;

INSTANTIATE_TEST_SUITE_P(IntegralsRectangleTests, DergynovSIntegralsFuncTest, kGtestValues, kFuncTestName);

}  // namespace
}  // namespace dergynov_s_integrals_multistep_rectangle
