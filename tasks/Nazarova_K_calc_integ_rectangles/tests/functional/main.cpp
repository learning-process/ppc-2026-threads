#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "Nazarova_K_calc_integ_rectangles/common/include/common.hpp"
#include "Nazarova_K_calc_integ_rectangles/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace nazarova_k_calc_integ_rectangles {

class NazarovaKCalcIntegRectanglesRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<3>(test_param);
  }

 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(param);
    expected_ = std::get<1>(param);
    eps_ = std::get<2>(param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_) <= eps_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_{0.0};
  double eps_{0.0};
};

namespace {

InType MakeInput(double (*function)(const std::vector<double> &), std::vector<double> lower_bounds,
                 std::vector<double> upper_bounds, std::vector<std::size_t> steps) {
  return InType{function, std::move(lower_bounds), std::move(upper_bounds), std::move(steps)};
}

double ConstFunction(const std::vector<double> & /*point*/) {
  return 5.0;
}

double Linear1D(const std::vector<double> &point) {
  return point[0];
}

double Product2D(const std::vector<double> &point) {
  return point[0] * point[1];
}

double Sum3D(const std::vector<double> &point) {
  return point[0] + point[1] + point[2];
}

double Square1D(const std::vector<double> &point) {
  return point[0] * point[0];
}

double Trig2D(const std::vector<double> &point) {
  return std::sin(point[0]) * std::cos(point[1]);
}

double ShiftedProduct3D(const std::vector<double> &point) {
  return (point[0] + 1.0) * (point[1] - 2.0) * point[2];
}

TEST_P(NazarovaKCalcIntegRectanglesRunFuncTests, MultidimensionalRectangleIntegration) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    TestType{MakeInput(ConstFunction, {0.0, 0.0}, {2.0, 3.0}, {8U, 9U}), 30.0, 1e-12, "Const2D"},
    TestType{MakeInput(Linear1D, {0.0}, {10.0}, {100U}), 50.0, 1e-12, "Linear1D"},
    TestType{MakeInput(Product2D, {0.0, 0.0}, {2.0, 3.0}, {20U, 30U}), 9.0, 1e-12, "Product2D"},
    TestType{MakeInput(Sum3D, {0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, {10U, 12U, 14U}), 18.0, 1e-12, "Sum3D"},
    TestType{MakeInput(Square1D, {0.0}, {1.0}, {1000U}), 1.0 / 3.0, 1e-7, "Square1D"},
    TestType{MakeInput(Trig2D, {0.0, 0.0}, {std::acos(-1.0), std::acos(-1.0) / 2.0}, {400U, 300U}), 2.0, 1e-5,
             "Trig2D"},
    TestType{MakeInput(ShiftedProduct3D, {-1.0, 2.0, 0.0}, {1.0, 4.0, 2.0}, {12U, 10U, 8U}), 8.0, 1e-12,
             "ShiftedProduct3D"}};

const auto kTestTasksList = ppc::util::AddFuncTask<NazarovaKCalcIntegRectanglesSEQ, InType>(
    kTestParam, PPC_SETTINGS_Nazarova_K_calc_integ_rectangles);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    NazarovaKCalcIntegRectanglesRunFuncTests::PrintFuncTestName<NazarovaKCalcIntegRectanglesRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RectangleIntegrationTests, NazarovaKCalcIntegRectanglesRunFuncTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace nazarova_k_calc_integ_rectangles
