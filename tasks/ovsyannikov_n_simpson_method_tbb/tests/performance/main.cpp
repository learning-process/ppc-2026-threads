#include <gtest/gtest.h>

#include <cmath>

#include "ovsyannikov_n_simpson_method_tbb/common/include/common.hpp"
#include "ovsyannikov_n_simpson_method_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace ovsyannikov_n_simpson_method_tbb {

class OvsyannikovNRunPerfTestTBB : public ppc::util::BaseRunPerfTests<InType, OutType> {
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
TEST_P(OvsyannikovNRunPerfTestTBB, SimpsonTestRunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, OvsyannikovNSimpsonMethodTBB>(PPC_SETTINGS_ovsyannikov_n_simpson_method_tbb);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = OvsyannikovNRunPerfTestTBB::CustomPerfTestName;
INSTANTIATE_TEST_SUITE_P(RunModeTests, OvsyannikovNRunPerfTestTBB, kGtestValues, kPerfTestName);
}  // namespace
}  // namespace ovsyannikov_n_simpson_method_tbb
