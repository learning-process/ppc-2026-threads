#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "shilin_n_monte_carlo_integration/common/include/common.hpp"
#include "shilin_n_monte_carlo_integration/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shilin_n_monte_carlo_integration {

class ShilinNRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &[lower, upper, n, func_type] = input_data_;
    double expected = IntegrandFunction::AnalyticalIntegral(func_type, lower, upper);

    double volume = 1.0;
    for (size_t i = 0; i < lower.size(); ++i) {
      volume *= (upper[i] - lower[i]);
    }
    double epsilon = std::max(volume * 10.0 / std::sqrt(static_cast<double>(n)), 1e-2);
    return std::abs(output_data - expected) <= epsilon;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(ShilinNRunFuncTestsThreads, MonteCarloIntegration) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {{
    std::make_tuple(std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, 10000, FuncType::kLinear),
                    "kLinear_1D"),
    std::make_tuple(std::make_tuple(std::vector<double>{0.0}, std::vector<double>{2.0}, 10000, FuncType::kSumSquares),
                    "kSumSquares_1D"),
    std::make_tuple(
        std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 50000, FuncType::kConstant),
        "kConstant_2D"),
    std::make_tuple(
        std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 50000, FuncType::kLinear),
        "kLinear_2D"),
    std::make_tuple(
        std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 50000, FuncType::kProduct),
        "kProduct_2D"),
    std::make_tuple(
        std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, 50000, FuncType::kSinProduct),
        "kSinProduct_2D"),
    std::make_tuple(std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, 100000,
                                    FuncType::kLinear),
                    "kLinear_3D"),
    std::make_tuple(std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, 100000,
                                    FuncType::kProduct),
                    "kProduct_3D"),
}};

const auto kTestTasksList = ppc::util::AddFuncTask<ShilinNMonteCarloIntegrationSEQ, InType>(
    kTestParam, PPC_SETTINGS_shilin_n_monte_carlo_integration);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ShilinNRunFuncTestsThreads::PrintFuncTestName<ShilinNRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MonteCarloTests, ShilinNRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shilin_n_monte_carlo_integration
