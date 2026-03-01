#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "../../common/include/common.hpp"
#include "../../seq/include/rect_method_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

class RectMethodFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string test_name = std::get<2>(params);

    std::cout << "[ DEBUG ] Test: " << test_name << '\n';
    std::cout << "[ DEBUG ] Expected: " << expected_output_ << '\n';
    std::cout << "[ DEBUG ] Actual:   " << output_data << '\n';
    std::cout << "[ DEBUG ] Diff:     " << std::abs(output_data - expected_output_) << '\n';
    return std::abs(output_data - expected_output_) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_{};
};

namespace {

TEST_P(RectMethodFuncTests, IntegrationTest) {
  ExecuteTest(GetParam());
}

IntegrationData test_1d = {
    .limits = {{0.0, 1.0}}, .n_steps = {100}, .func = [](const std::vector<double> &x) { return x[0]; }};

IntegrationData test_2d = {.limits = {{0.0, 1.0}, {0.0, 1.0}},
                           .n_steps = {50, 50},
                           .func = [](const std::vector<double> &x) { return x[0] + x[1]; }};

IntegrationData test_3d = {.limits = {{0.0, 2.0}, {0.0, 3.0}, {0.0, 1.0}},
                           .n_steps = {20, 20, 20},
                           .func = [](const std::vector<double> & /*unused*/) { return 1.0; }};

const std::array<TestType, 3> kTestCases = {std::make_tuple(test_1d, 0.5, "1D_Linear"),
                                            std::make_tuple(test_2d, 1.0, "2D_Linear"),
                                            std::make_tuple(test_3d, 6.0, "3D_Constant")};

const auto kTestTasksList = ppc::util::AddFuncTask<RectMethodSequential, InType>(
    kTestCases, PPC_SETTINGS_kutergin_v_multidimensional_integration_rect_method);

const auto kGTestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = RectMethodFuncTests::PrintFuncTestName<RectMethodFuncTests>;

INSTANTIATE_TEST_SUITE_P(MultidimensionalIntegration, RectMethodFuncTests, kGTestValues, kTestName);

}  // namespace

}  // namespace kutergin_v_multidimensional_integration_rect_method
