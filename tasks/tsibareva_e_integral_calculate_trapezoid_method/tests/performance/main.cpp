#include <gtest/gtest.h>

//#include "tsibareva_e_integral_calculate_trapezoid_method/all/include/ops_all.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/common/include/common.hpp"
//#include "tsibareva_e_integral_calculate_trapezoid_method/omp/include/ops_omp.hpp"
#include "tsibareva_e_integral_calculate_trapezoid_method/seq/include/ops_seq.hpp"
//#include "tsibareva_e_integral_calculate_trapezoid_method/stl/include/ops_stl.hpp"
//#include "tsibareva_e_integral_calculate_trapezoid_method/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

class TsibarevaERunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    input_data_.dimension = 3;
    input_data_.lower_bounds = {0.0, 0.0, 0.0};
    input_data_.upper_bounds = {1.0, 1.0, 1.0};
    input_data_.num_steps = {200, 200, 200};
    input_data_.function = [](const std::vector<double>& x) { 
      return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; 
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    double expected = 1.0;
    return std::fabs(output_data - expected) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TsibarevaERunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, 
    /*TsibarevaEIntegralCalculateTrapezoidMethodALL, */
    /*TsibarevaEIntegralCalculateTrapezoidMethodOMP, */
    TsibarevaEIntegralCalculateTrapezoidMethodSEQ>      
    /*TsibarevaEIntegralCalculateTrapezoidMethodSTL,*/ 
    /*TsibarevaEIntegralCalculateTrapezoidMethodTBB>,*/
    (PPC_SETTINGS_tsibareva_e_integral_calculate_trapezoid_method);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TsibarevaERunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TsibarevaERunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsibareva_e_integral_calculate_trapezoid_method