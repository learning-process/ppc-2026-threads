#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/common/include/common.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/omp/include/ops_omp.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/seq/include/ops_seq.hpp"
#include "vinyaikina_e_multidimensional_integrals_simpson_method/tbb/include/ops_tbb.hpp"

namespace vinyaikina_e_multidimensional_integrals_simpson_method {
namespace {

class VinyaikinaESimpsonFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_ = std::get<1>(params);
    etalon_ = std::get<2>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-3;
    return std::fabs(output_data - etalon_) <= eps;
  }

  InType GetTestInputData() final {
    return input_;
  }

 private:
  InType input_;
  OutType etalon_ = 0.0;
};

TEST_P(VinyaikinaESimpsonFuncTests, Run) {
  ExecuteTest(GetParam());
}

double CountNDimArea(const std::vector<std::pair<double, double>> &borders) {
  double area = 1.0;
  for (const auto &border : borders) {
    area *= (border.second - border.first);
  }
  return area;
}

double IntLinear1d(double a, double b) {
  return (b * b - a * a) / 2.0;
}

double Intsin1d(double a, double b) {
  return std::cos(a) - std::cos(b);
}

double Intx21d(double a, double b) {
  return (b * b * b - a * a * a) / 3.0;
}

double Intx31d(double a, double b) {
  return (b * b * b * b - a * a * a * a) / 4.0;
}

double Intx41d(double a, double b) {
  return (b * b * b * b * b - a * a * a * a * a) / 5.0;
}

double IntExp1d(double a, double b) {
  return std::exp(b) - std::exp(a);
}

double Intcos1d(double a, double b) {
  return std::sin(b) - std::sin(a);
}

double Intxy2d(double a1, double b1, double a2, double b2) {
  return ((b1 * b1 - a1 * a1) * (b2 * b2 - a2 * a2)) / 4.0;
}

double IntxPlusy2d(double a1, double b1, double a2, double b2) {
  return ((b2 - a2) * (b1 * b1 - a1 * a1) / 2.0) + ((b1 - a1) * (b2 * b2 - a2 * a2) / 2.0);
}

double Intx2y2d(double a1, double b1, double a2, double b2) {
  return Intx21d(a1, b1) * IntLinear1d(a2, b2);
}

double Intxyz3d(double a1, double b1, double a2, double b2, double a3, double b3) {
  return IntLinear1d(a1, b1) * IntLinear1d(a2, b2) * IntLinear1d(a3, b3);
}

double Intx2Y2Z23d(double a1, double b1, double a2, double b2, double a3, double b3) {
  return Intx21d(a1, b1) * Intx21d(a2, b2) * Intx21d(a3, b3);
}

double IntexpSum3d(double a1, double b1, double a2, double b2, double a3, double b3) {
  return IntExp1d(a1, b1) * IntExp1d(a2, b2) * IntExp1d(a3, b3);
}

auto one = [](const std::vector<double> &) { return 1.0; };
auto linear1d = [](const std::vector<double> &x) { return x[0]; };
auto sin1d = [](const std::vector<double> &x) { return std::sin(x[0]); };
auto x21d = [](const std::vector<double> &x) { return x[0] * x[0]; };
auto x31d = [](const std::vector<double> &x) { return x[0] * x[0] * x[0]; };
auto x41d = [](const std::vector<double> &x) { return std::pow(x[0], 4); };
auto exp1d = [](const std::vector<double> &x) { return std::exp(x[0]); };
auto cos1d = [](const std::vector<double> &x) { return std::cos(x[0]); };

auto xy2d = [](const std::vector<double> &x) { return x[0] * x[1]; };
auto x_plus_y2d = [](const std::vector<double> &x) { return x[0] + x[1]; };
auto x2_y2d = [](const std::vector<double> &x) { return x[0] * x[0] * x[1]; };

auto xyz_3d = [](const std::vector<double> &x) { return x[0] * x[1] * x[2]; };
auto x2_y2_z2_3d = [](const std::vector<double> &x) { return x[0] * x[0] * x[1] * x[1] * x[2] * x[2]; };
auto exp_sum_3d = [](const std::vector<double> &x) { return std::exp(x[0] + x[1] + x[2]); };

const std::array<TestType, 16> kTests = {{
    TestType{"area1d_0_1", InType{0.005, {{0.0, 1.0}}, one}, CountNDimArea({{0.0, 1.0}})},

    TestType{"area2d_0_1_x_0_1", InType{0.005, {{0.0, 1.0}, {0.0, 1.0}}, one}, CountNDimArea({{0.0, 1.0}, {0.0, 1.0}})},

    TestType{"volume_3d_0_05_x3", InType{0.005, {{0.0, 0.25}, {0.0, 0.25}, {0.0, 0.25}}, one},
             CountNDimArea({{0.0, 0.25}, {0.0, 0.25}, {0.0, 0.25}})},

    TestType{"Linear1d_0_2", InType{0.01, {{0.0, 2.0}}, linear1d}, IntLinear1d(0.0, 2.0)},

    TestType{"x21d_0_1", InType{0.005, {{0.0, 1.0}}, x21d}, Intx21d(0.0, 1.0)},

    TestType{"x31d_0_1", InType{0.005, {{0.0, 1.0}}, x31d}, Intx31d(0.0, 1.0)},

    TestType{"x41d_0_1", InType{0.005, {{0.0, 1.0}}, x41d}, Intx41d(0.0, 1.0)},

    TestType{"exp1d_0_1", InType{0.005, {{0.0, 1.0}}, exp1d}, IntExp1d(0.0, 1.0)},

    TestType{"cos1d_0_pi2", InType{0.001, {{0.0, std::numbers::pi / 2.0}}, cos1d},
             Intcos1d(0.0, std::numbers::pi / 2.0)},

    TestType{"sin1d_0_pi", InType{0.001, {{0.0, std::numbers::pi}}, sin1d}, Intsin1d(0.0, std::numbers::pi)},

    TestType{"xy2d_0_1_x_0_1", InType{0.005, {{0.0, 1.0}, {0.0, 1.0}}, xy2d}, Intxy2d(0.0, 1.0, 0.0, 1.0)},

    TestType{"x_plus_y2d_0_1_x_0_1", InType{0.005, {{0.0, 1.0}, {0.0, 1.0}}, x_plus_y2d},
             IntxPlusy2d(0.0, 1.0, 0.0, 1.0)},

    TestType{"x2_y2d_0_1_x_0_1", InType{0.005, {{0.0, 1.0}, {0.0, 1.0}}, x2_y2d}, Intx2y2d(0.0, 1.0, 0.0, 1.0)},

    TestType{"xyz_3d_0_1_x3", InType{0.005, {{0.0, 0.25}, {0.0, 0.25}, {0.0, 0.25}}, xyz_3d},
             Intxyz3d(0.0, 0.25, 0.0, 0.25, 0.0, 0.25)},

    TestType{"x2_y2_z2_3d_0_1_x3", InType{0.01, {{0.0, 0.25}, {0.0, 0.25}, {0.0, 0.25}}, x2_y2_z2_3d},
             Intx2Y2Z23d(0.0, 0.25, 0.0, 0.25, 0.0, 0.25)},

    TestType{"exp_sum_3d_0_05_x3", InType{0.005, {{0.0, 0.25}, {0.0, 0.25}, {0.0, 0.25}}, exp_sum_3d},
             IntexpSum3d(0.0, 0.25, 0.0, 0.25, 0.0, 0.25)},
}};

const auto kTaskName = PPC_SETTINGS_vinyaikina_e_multidimensional_integrals_simpson_method;

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<VinyaikinaEMultidimIntegrSimpsonSEQ, InType>(kTests, kTaskName),
                   ppc::util::AddFuncTask<VinyaikinaEMultidimIntegrSimpsonOMP, InType>(kTests, kTaskName),
                   ppc::util::AddFuncTask<VinyaikinaEMultidimIntegrSimpsonTBB, InType>(kTests, kTaskName));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = VinyaikinaESimpsonFuncTests::PrintFuncTestName<VinyaikinaESimpsonFuncTests>;
;

INSTANTIATE_TEST_SUITE_P(MultidinIntegralsSimpsonTests, VinyaikinaESimpsonFuncTests, kGtestValues, kFuncTestName);

}  // namespace
}  // namespace vinyaikina_e_multidimensional_integrals_simpson_method
