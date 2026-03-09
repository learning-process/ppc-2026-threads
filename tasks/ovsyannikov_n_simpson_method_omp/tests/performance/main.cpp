#include <gtest/gtest.h>

#include <cmath>

#include "ovsyannikov_n_simpson_method_omp/common/include/common.hpp"
#include "ovsyannikov_n_simpson_method_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace ovsyannikov_n_simpson_method_omp {

class OvsyannikovNRunPerfTestThreadsOMP : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    input_data_ = InType{0.0, 1.0, 0.0, 1.0, 2000, 2000};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - 1.0) < 1e-4;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {};
};

namespace {
TEST_P(OvsyannikovNRunPerfTestThreadsOMP, SimpsonTestRunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, OvsyannikovNSimpsonMethodOMP>(PPC_SETTINGS_ovsyannikov_n_simpson_method_omp);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = OvsyannikovNRunPerfTestThreadsOMP::CustomPerfTestName;
INSTANTIATE_TEST_SUITE_P(RunModeTests, OvsyannikovNRunPerfTestThreadsOMP, kGtestValues, kPerfTestName);
}  // namespace
}  // namespace ovsyannikov_n_simpson_method_omp
