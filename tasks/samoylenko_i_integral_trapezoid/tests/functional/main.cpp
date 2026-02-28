#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "samoylenko_i_integral_trapezoid/common/include/common.hpp"
#include "samoylenko_i_integral_trapezoid/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace samoylenko_i_integral_trapezoid {

class SamoylenkoIIntegralTrapezoidFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &in = std::get<0>(test_param);
    return std::to_string(in.a.size()) + "_dims__" + std::to_string(in.function_choice) + "_func";
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) override {
    return std::abs(output_data - expected_) < 1e-3;
  }

  InType GetTestInputData() override {
    return input_data_;
  }

 private:
  InType input_data_;
  double expected_ = 0.0;
};

namespace {

TEST_P(SamoylenkoIIntegralTrapezoidFuncTests, IntegrationTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_pair(InType{{0.0}, {10.0}, {1000}, 0}, 50.0),
    std::make_pair(InType{{0.0, 0.0}, {2.0, 2.0}, {100, 100}, 0}, 8.0),
    std::make_pair(InType{{0.0, 0.0}, {2.0, 3.0}, {100, 100}, 1}, 9.0),
    std::make_pair(InType{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {50, 50, 50}, 2}, 1.0),
    std::make_pair(InType{{0.0}, {std::numbers::pi}, {1000}, 3}, 2.0)};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<SamoylenkoIIntegralTrapezoidSEQ, InType>(
    kTestParam, PPC_SETTINGS_samoylenko_i_integral_trapezoid));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName =
    SamoylenkoIIntegralTrapezoidFuncTests::PrintFuncTestName<SamoylenkoIIntegralTrapezoidFuncTests>;

INSTANTIATE_TEST_SUITE_P(BasicTests, SamoylenkoIIntegralTrapezoidFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace samoylenko_i_integral_trapezoid
