#include <gtest/gtest.h>

#include "ovsyannikov_n_simpson_method/common/include/common.hpp"
#include "ovsyannikov_n_simpson_method/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace ovsyannikov_n_simpson_method {

class OvsyannikovNRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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
  InType input_data_;
};

TEST_P(OvsyannikovNRunPerfTestThreads, SimpsonTestRunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, OvsyannikovNSimpsonMethodSEQ>(PPC_SETTINGS_ovsyannikov_n_simpson_method);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = OvsyannikovNRunPerfTestThreads::CustomPerfTestName;
INSTANTIATE_TEST_SUITE_P(RunModeTests, OvsyannikovNRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace ovsyannikov_n_simpson_method
