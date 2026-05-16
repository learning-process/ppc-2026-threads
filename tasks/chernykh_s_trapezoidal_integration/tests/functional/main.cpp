#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"
#include "chernykh_s_trapezoidal_integration/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace chernykh_s_trapezoidal_integration {

class ChernykhSRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_output_) < 1e-3;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{{}, {}, nullptr};
  OutType expected_output_ = 0.0;
};

namespace {

TEST_P(ChernykhSRunFuncTestsThreads, SEQIntegration) {
  ExecuteTest(GetParam());
}

double FConst(const std::vector<double> & /*unused*/) {
  return 5.0;
}
double FLinear(const std::vector<double> &x) {
  return x[0];
}
double FSum2d(const std::vector<double> &x) {
  return x[0] + x[1];
}
double FParabola(const std::vector<double> &x) {
  return x[0] * x[0];
}
double FSin(const std::vector<double> &x) {
  return std::sin(x[0]);
}

const std::array<TestType, 6> kTestParam = {
    std::make_tuple(InType({{0.0, 1.0}}, {1000}, FLinear), 0.5, "Linear_1D"),
    std::make_tuple(InType({{0.0, 1.0}, {0.0, 1.0}}, {100, 100}, FSum2d), 1.0, "Sum_2D"),
    std::make_tuple(InType({{0.0, 2.0}, {0.0, 2.0}, {0.0, 2.0}}, {10, 10, 10}, FConst), 40.0, "Constant_3D"),
    std::make_tuple(InType({{0.0, 1.0}}, {1000}, FParabola), 0.3333333333, "Parabola_1D"),
    std::make_tuple(InType({{0.0, std::numbers::pi}}, {1000}, FSin), 2.0, "Sin_1D"),
    std::make_tuple(InType({{1.0, 1.0}}, {100}, FLinear), 0.0, "Zero_Range")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<ChernykhSTrapezoidalIntegrationSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ChernykhSRunFuncTestsThreads::PrintFuncTestName<ChernykhSRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(ChernykhSIntegrationTests, ChernykhSRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace chernykh_s_trapezoidal_integration
