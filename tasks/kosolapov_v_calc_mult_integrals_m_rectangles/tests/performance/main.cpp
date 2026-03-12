#include <gtest/gtest.h>

#include <cmath>
#include <tuple>

#include "kosolapov_v_calc_mult_integrals_m_rectangles/common/include/common.hpp"
#include "kosolapov_v_calc_mult_integrals_m_rectangles/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kosolapov_v_calc_mult_integrals_m_rectangles {

class KosolapovVCalcMultIntegralsMRectanglesPerfTestProcesses : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 10000;
  InType input_data_;

  void SetUp() override {
    int func_id = 1;
    input_data_ = std::make_tuple(kCount_, func_id);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - (2.0 / 3.0)) < 0.01;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KosolapovVCalcMultIntegralsMRectanglesPerfTestProcesses, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KosolapovVCalcMultIntegralsMRectanglesSEQ>(
    PPC_SETTINGS_kosolapov_v_calc_mult_integrals_m_rectangles);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KosolapovVCalcMultIntegralsMRectanglesPerfTestProcesses::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunRectanglesModeTests, KosolapovVCalcMultIntegralsMRectanglesPerfTestProcesses, kGtestValues,
                         kPerfTestName);

}  // namespace kosolapov_v_calc_mult_integrals_m_rectangles
