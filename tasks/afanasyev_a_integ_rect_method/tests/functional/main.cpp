#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"
#include "afanasyev_a_integ_rect_method/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace afanasyev_a_integ_rect_method {

class AfanasyevARunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    n_ = std::get<0>(params);
    if (n_ >= 80) {
      tol_ = 2e-4;
    } else if (n_ >= 30) {
      tol_ = 1e-3;
    } else {
      tol_ = 1e-2;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    constexpr int kDim = 3;
    const double pi = std::acos(-1.0);
    const double i1 = 0.5 * std::sqrt(pi) * std::erf(1.0);
    const double expected = std::pow(i1, kDim);
    const double err = std::fabs(output_data - expected);
    return err <= tol_;
  }

  InType GetTestInputData() final {
    return n_;
  }

 private:
  InType n_ = 10;
  double tol_ = 1e-2;
};

namespace {

TEST_P(AfanasyevARunFuncTests, MidpointApproximation) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(10, "n10"), std::make_tuple(30, "n30"),
                                            std::make_tuple(80, "n80")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<AfanasyevAIntegRectMethodSEQ, InType>(
    kTestParam, PPC_SETTINGS_afanasyev_a_integ_rect_method));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = AfanasyevARunFuncTests::PrintFuncTestName<AfanasyevARunFuncTests>;

INSTANTIATE_TEST_SUITE_P(MidpointTests, AfanasyevARunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace afanasyev_a_integ_rect_method
