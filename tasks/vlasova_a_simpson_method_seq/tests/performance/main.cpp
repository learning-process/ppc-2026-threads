#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "vlasova_a_simpson_method_seq/common/include/common.hpp"
#include "vlasova_a_simpson_method_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace vlasova_a_simpson_method_seq {

double Gaussian3D(const std::vector<double>& x) {
  return std::exp(-(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
}

class VlasovaASimpsonMethodPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    std::vector<double> a = {-2.0, -2.0, -2.0};
    std::vector<double> b = {2.0, 2.0, 2.0};
    std::vector<int> n = {200, 200, 200};
    
    input_data_ = SimpsonTask(Gaussian3D, a, b, n);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::isfinite(output_data) && output_data > 0.0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {
TEST_P(VlasovaASimpsonMethodPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, VlasovaASimpsonMethodSEQ>(PPC_SETTINGS_vlasova_a_simpson_method_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = VlasovaASimpsonMethodPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SimpsonMethodPerfTests, VlasovaASimpsonMethodPerfTests, kGtestValues, kPerfTestName);

} // namespace

}  // namespace vlasova_a_simpson_method_seq