#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "shkrebko_m_calc_of_integral_rect/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

class ShkrebkoMRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(ShkrebkoMRunFuncTests, MultiDimRectangleMethod) {
  ExecuteTest(GetParam());
}
const std::array<TestType, 5> kTestParam = {
    TestType{"Const_1D_0_1_N100", InType{{{0.0, 1.0}}, 100, [](const std::vector<double> &) { return 1.0; }}, 1.0},
    TestType{"Linear_1D_0_2_N200", InType{{{0.0, 2.0}}, 200, [](const std::vector<double> &x) { return x[0]; }}, 2.0},
    TestType{"Quad_1D_m1_1_N150", InType{{{-1.0, 1.0}}, 150, [](const std::vector<double> &x) { return x[0] * x[0]; }},
             2.0 / 3.0},
    TestType{"Prod_2D_0_2_1_3_N100",
             InType{{{0.0, 2.0}, {1.0, 3.0}}, 100, [](const std::vector<double> &x) { return x[0] * x[1]; }}, 8.0},
    TestType{"Trig_2D_sin_sin_pi_pi2_N200",
             InType{{{0.0, std::numbers::pi}, {0.0, std::numbers::pi / 2.0}},
                    200,
                    [](const std::vector<double> &x) { return std::sin(x[0]) * std::sin(x[1]); }},
             2.0},
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ShkrebkoMCalcOfIntegralRectSEQ, InType>(
    kTestParam, PPC_SETTINGS_shkrebko_m_calc_of_integral_rect));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = ShkrebkoMRunFuncTests::PrintFuncTestName<ShkrebkoMRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(IntegralsRectangleMethodTests, ShkrebkoMRunFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace shkrebko_m_calc_of_integral_rect
