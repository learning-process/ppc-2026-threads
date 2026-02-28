#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaERunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    auto params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    IntegralTestType test_type = std::get<0>(params);
    
    input_data_ = GenerateIntegralInput(test_type);
    expected_output_ = GenerateExpectedOutput(test_type);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    double tolerance = 1e-4;
    return std::fabs(output_data - expected_output_) < tolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  double expected_output_;
};

namespace {

TEST_P(TsibarevaERunFuncTestsThreads, IntegralCalculation) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParams = {
    std::make_tuple(IntegralTestType::SUCCESS_SIMPLE_2D, "2d_simple"),
    std::make_tuple(IntegralTestType::SUCCESS_CONSTANT_2D, "2d_constant"),
    std::make_tuple(IntegralTestType::SUCCESS_SIMPLE_3D, "3d_simple"),
    std::make_tuple(IntegralTestType::SUCCESS_CONSTANT_3D, "3d_constant"),
    std::make_tuple(IntegralTestType::INVALID_LOWER_BOUND_EQUAL, "invalid_lower_bound_equal"),
    std::make_tuple(IntegralTestType::INVALID_STEPS_NEGATIVE, "invalid_steps_negative"),
    std::make_tuple(IntegralTestType::INVALID_EMPTY_BOUNDS, "invalid_empty_bounds")
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TsibarevaEIntegralCalculateTrapezoidMethodSEQ, InType>(kTestParams, PPC_SETTINGS_tsibareva_e_integral_calculate_trapezoid_method)
);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = TsibarevaERunFuncTestsThreads::PrintFuncTestName<TsibarevaERunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(IntegralCalculation, TsibarevaERunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsibareva_e_integral_calculate_trapezoid_method