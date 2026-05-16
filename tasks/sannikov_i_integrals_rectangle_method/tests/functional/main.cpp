#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"
#include "sannikov_i_integrals_rectangle_method/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sannikov_i_integrals_rectangle_method {
class SannikovIRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
    const double eps = 1e-4;
    return std::fabs(output_data - expected_) <= eps;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_ = 0.0;
};

namespace {

TEST_P(SannikovIRunFuncTests, MultiDimRectangleMethod) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 23> kTestParam = {
    TestType{"Const_1D_0_1_N100", InType{[](const std::vector<double> &) { return 1.0; }, {{0.0, 1.0}}, 100}, 1.0},
    TestType{"Const_1D_m2_3_N200", InType{[](const std::vector<double> &) { return 1.0; }, {{-2.0, 3.0}}, 200}, 5.0},
    TestType{"Const_2D_0_1_x_0_2_N80",
             InType{[](const std::vector<double> &) { return 1.0; }, {{0.0, 1.0}, {0.0, 2.0}}, 80}, 2.0},
    TestType{"Const_3D_0_1_x_m1_1_x_2_5_N40",
             InType{[](const std::vector<double> &) { return 1.0; }, {{0.0, 1.0}, {-1.0, 1.0}, {2.0, 5.0}}, 40}, 6.0},
    TestType{"Const_4D_0_2_x_m1_1_x_2_3_x_0_05_N30", InType{[](const std::vector<double> &) {
  return 1.0;
}, {{0.0, 2.0}, {-1.0, 1.0}, {2.0, 3.0}, {0.0, 0.5}}, 30},
             2.0},
    TestType{"Linear_1D_x_0_2_N200", InType{[](const std::vector<double> &x) { return x[0]; }, {{0.0, 2.0}}, 200}, 2.0},
    TestType{"Linear_1D_x_m1_1_N200", InType{[](const std::vector<double> &x) { return x[0]; }, {{-1.0, 1.0}}, 200},
             0.0},
    TestType{"Linear_2D_x_plus_y_0_1_x_0_1_N80",
             InType{[](const std::vector<double> &x) { return x[0] + x[1]; }, {{0.0, 1.0}, {0.0, 1.0}}, 80}, 1.0},
    TestType{"Linear_3D_sum_0_1_x_0_2_x_m1_1_N80", InType{[](const std::vector<double> &x) {
  return x[0] + x[1] + x[2];
}, {{0.0, 1.0}, {0.0, 2.0}, {-1.0, 1.0}}, 80},
             6.0},
    TestType{"Quad_1D_x2_m1_2_N200",
             InType{[](const std::vector<double> &x) { return x[0] * x[0]; }, {{-1.0, 2.0}}, 200}, 3.0},
    TestType{"Quad_2D_x2_plus_y2_0_1_x_0_2_N100", InType{[](const std::vector<double> &x) {
  return (x[0] * x[0]) + (x[1] * x[1]);
}, {{0.0, 1.0}, {0.0, 2.0}}, 100},
             10.0 / 3.0},
    TestType{"Quad_1D_x2_m1_1_N150",
             InType{[](const std::vector<double> &x) { return x[0] * x[0]; }, {{-1.0, 1.0}}, 150}, 2.0 / 3.0},
    TestType{"Prod_2D_xy_0_2_x_1_3_N100",
             InType{[](const std::vector<double> &x) { return x[0] * x[1]; }, {{0.0, 2.0}, {1.0, 3.0}}, 100}, 8.0},
    TestType{"Prod_3D_xyz_0_1_x_0_2_x_0_3_N100", InType{[](const std::vector<double> &x) {
  return x[0] * x[1] * x[2];
}, {{0.0, 1.0}, {0.0, 2.0}, {0.0, 3.0}}, 100},
             4.5},
    TestType{"Prod_2D_xy_m1_1_x_0_2_N100",
             InType{[](const std::vector<double> &x) { return x[0] * x[1]; }, {{-1.0, 1.0}, {0.0, 2.0}}, 100}, 0.0},
    TestType{"Trig_1D_sin_0_pi_N200",
             InType{[](const std::vector<double> &x) { return std::sin(x[0]); }, {{0.0, std::numbers::pi}}, 200}, 2.0},
    TestType{"Trig_1D_cos_0_pi_N200",
             InType{[](const std::vector<double> &x) { return std::cos(x[0]); }, {{0.0, std::numbers::pi}}, 200}, 0.0},
    TestType{"Trig_2D_sin_sin_0_pi_x_0_pi2_N200", InType{[](const std::vector<double> &x) {
  return std::sin(x[0]) * std::sin(x[1]);
}, {{0.0, std::numbers::pi}, {0.0, std::numbers::pi / 2.0}}, 200},
             2.0},
    TestType{"Trig_1D_cos_0_pi2_N200",
             InType{[](const std::vector<double> &x) { return std::cos(x[0]); }, {{0.0, std::numbers::pi / 2.0}}, 200},
             1.0},
    TestType{"Exp_1D_exp_0_1_N200",
             InType{[](const std::vector<double> &x) { return std::exp(x[0]); }, {{0.0, 1.0}}, 200},
             std::numbers::e - 1.0},
    TestType{"Exp_2D_exp_x_plus_y_0_1_x_0_1_N200",
             InType{[](const std::vector<double> &x) { return std::exp(x[0] + x[1]); }, {{0.0, 1.0}, {0.0, 1.0}}, 200},
             (std::numbers::e - 1.0) * (std::numbers::e - 1.0)},
    TestType{"Const_1D_small_0_001_N100", InType{[](const std::vector<double> &) { return 1.0; }, {{0.0, 1e-3}}, 100},
             1e-3},
    TestType{"Const_2D_scale_0_100_x_0_001_N200",
             InType{[](const std::vector<double> &) { return 1.0; }, {{0.0, 100.0}, {0.0, 0.01}}, 200}, 1.0},
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<SannikovIIntegralsRectangleMethodSEQ, InType>(
    kTestParam, PPC_SETTINGS_sannikov_i_integrals_rectangle_method));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = SannikovIRunFuncTests::PrintFuncTestName<SannikovIRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(IntegralsRectangleMethodTests, SannikovIRunFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace sannikov_i_integrals_rectangle_method
