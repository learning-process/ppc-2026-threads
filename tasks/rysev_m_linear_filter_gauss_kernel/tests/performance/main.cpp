#include <gtest/gtest.h>

#include "rysev_m_linear_filter_gauss_kernel/common/include/common.hpp"
#include "rysev_m_linear_filter_gauss_kernel/omp/include/ops_omp.hpp"
#include "rysev_m_linear_filter_gauss_kernel/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace rysev_m_linear_filter_gauss_kernel {

class RysevMPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kImageSize_ = 512;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kImageSize_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data > 0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(RysevMPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RysevMGaussFilterSEQ, RysevMGaussFilterOMP>(
    PPC_SETTINGS_rysev_m_linear_filter_gauss_kernel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = RysevMPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RysevMPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace rysev_m_linear_filter_gauss_kernel
