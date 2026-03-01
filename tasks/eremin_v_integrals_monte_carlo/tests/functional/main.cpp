#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"
#include "eremin_v_integrals_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace eremin_v_integrals_monte_carlo {

class EreminVRunFuncTestsThreadsIntegralsMonteCarlo : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &input = std::get<0>(test_param);
    std::string result = "dims_" + std::to_string(input.bounds.size()) + "_samples_" + std::to_string(input.samples);
    std::ranges::replace(result, '.', '_');
    std::ranges::replace(result, '-', 'm');
    return result;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = std::get<0>(params);

    expected_result_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    double tolerance = 1e-2;
    return std::abs(output_data - expected_result_) <= tolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_{};
};

namespace {

TEST_P(EreminVRunFuncTestsThreadsIntegralsMonteCarlo, IntegralsMonteCarloFunc) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {
    std::make_tuple(MonteCarloInput{.bounds = {{0.0, 1.0}},  // 1D
                                    .samples = 1'000'000,
                                    .func = [](const std::vector<double> &x) { return x[0] * x[0]; }},
                    1.0 / 3.0),

    std::make_tuple(MonteCarloInput{.bounds = {{0.0, 1.0}, {0.0, 1.0}},  // 2D
                                    .samples = 1'000'000,
                                    .func = [](const std::vector<double> &x) { return x[0] * x[1]; }},
                    0.25),

    std::make_tuple(MonteCarloInput{.bounds = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}},  // 3D
                                    .samples = 1'000'000,
                                    .func = [](const std::vector<double> &x) { return x[0] + x[1] + x[2]; }},
                    1.5)};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<EreminVIntegralsMonteCarloSEQ, InType>(
    kTestParam, PPC_SETTINGS_eremin_v_integrals_monte_carlo));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    EreminVRunFuncTestsThreadsIntegralsMonteCarlo::PrintFuncTestName<EreminVRunFuncTestsThreadsIntegralsMonteCarlo>;

INSTANTIATE_TEST_SUITE_P(IntegralsMonteCarloTestsFunc, EreminVRunFuncTestsThreadsIntegralsMonteCarlo, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace eremin_v_integrals_monte_carlo
