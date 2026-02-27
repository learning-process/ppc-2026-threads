#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "shkrebko_m_calc_of_integral_rect/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

class ShkrebkoMRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_ = std::get<1>(params);
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

InType MakeInType(std::function<double(const std::vector<double> &)> func,
                  const std::vector<std::pair<double, double>> &limits, int n_steps) {
  std::vector<int> steps(limits.size(), n_steps);
  return InType{limits, steps, func};
}

const std::array<TestType, 10> kTestCases = {
    TestType{MakeInType([](const std::vector<double> &) { return 1.0; }, {{0.0, 1.0}}, 100), 1.0, "Const_1D"},
    TestType{MakeInType([](const std::vector<double> &) { return 1.0; }, {{0.0, 1.0}, {0.0, 2.0}}, 80), 2.0,
             "Const_2D"},
    TestType{MakeInType([](const std::vector<double> &) { return 1.0; }, {{0.0, 1.0}, {-1.0, 1.0}, {2.0, 5.0}}, 40),
             6.0, "Const_3D"},
    TestType{MakeInType([](const std::vector<double> &x) { return x[0]; }, {{0.0, 2.0}}, 200), 2.0, "Linear_1D"},
    TestType{MakeInType([](const std::vector<double> &x) { return x[0] + x[1]; }, {{0.0, 1.0}, {0.0, 1.0}}, 80), 1.0,
             "Linear_2D"},
    TestType{MakeInType([](const std::vector<double> &x) { return x[0] * x[0]; }, {{-1.0, 2.0}}, 200), 3.0, "Quad_1D"},
    TestType{MakeInType([](const std::vector<double> &x) { return x[0] * x[1]; }, {{0.0, 2.0}, {1.0, 3.0}}, 100), 8.0,
             "Prod_2D"},
    TestType{MakeInType([](const std::vector<double> &x) { return std::sin(x[0]); }, {{0.0, std::numbers::pi}}, 200),
             2.0, "Trig_sin_1D"},
    TestType{MakeInType([](const std::vector<double> &x) { return std::exp(x[0]); }, {{0.0, 1.0}}, 200),
             std::numbers::e - 1.0, "Exp_1D"},
    TestType{MakeInType([](const std::vector<double> &) { return 1.0; }, {{0.0, 1e-3}}, 100), 1e-3, "Const_1D_small"},
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ShkrebkoMCalcOfIntegralRectSEQ, InType>(
    kTestCases, PPC_SETTINGS_shkrebko_m_calc_of_integral_rect));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = ShkrebkoMRunFuncTests::PrintFuncTestName<ShkrebkoMRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(IntegralsRectangleMethodTests, ShkrebkoMRunFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace shkrebko_m_calc_of_integral_rect
