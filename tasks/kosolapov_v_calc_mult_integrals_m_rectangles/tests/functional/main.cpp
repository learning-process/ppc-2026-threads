#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>

#include "kosolapov_v_calc_mult_integrals_m_rectangles/common/include/common.hpp"
#include "kosolapov_v_calc_mult_integrals_m_rectangles/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kosolapov_v_calc_mult_integrals_m_rectangles {

class KosolapovVCalcMultIntegralsMRectanglesFuncTestsProcesses
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_count_steps_" + std::to_string(std::get<1>(test_param)) +
           "_func_id_" + std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int steps = std::get<0>(params);
    int func_id = std::get<1>(params);
    input_data_ = std::make_tuple(steps, func_id);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int func_id = std::get<1>(input_data_);
    switch (func_id) {
      case 1:
        return std::abs(output_data - (2.0 / 3.0)) < 0.01;
      case 2:
        return std::abs(output_data - 2.0) < 0.05;
      case 3:
        return std::abs(output_data - 2.230) < 0.05;
      case 4:
        return std::abs(output_data) < 0.05;
      default:
        return false;
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KosolapovVCalcMultIntegralsMRectanglesFuncTestsProcesses, CalcMultIntegralsMRectangles) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {std::make_tuple(20, 1, "test"), std::make_tuple(20, 2, "test"),
                                            std::make_tuple(20, 3, "test"), std::make_tuple(20, 4, "test"),
                                            std::make_tuple(40, 1, "test"), std::make_tuple(40, 2, "test"),
                                            std::make_tuple(40, 3, "test"), std::make_tuple(40, 4, "test")};

const auto kTestTasksList = std::tuple_cat(
                                           ppc::util::AddFuncTask<KosolapovVCalcMultIntegralsMRectanglesSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_kosolapov_v_calc_mult_integrals_m_rectangles));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KosolapovVCalcMultIntegralsMRectanglesFuncTestsProcesses::PrintFuncTestName<
    KosolapovVCalcMultIntegralsMRectanglesFuncTestsProcesses>;

INSTANTIATE_TEST_SUITE_P(CalcMultIntegralsMRectangles, KosolapovVCalcMultIntegralsMRectanglesFuncTestsProcesses, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace kosolapov_v_calc_mult_integrals_m_rectangles
