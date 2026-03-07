#include <gtest/gtest.h>

#include "gasenin_l_djstra_omp/common/include/common.hpp"
#include "gasenin_l_djstra_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraOmpPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};
  OutType expected_output_{};

  void SetUp() override {
    input_data_ = 200;
    expected_output_ = input_data_ * (input_data_ - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(GaseninLDjstraOmpPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GaseninLDjstraOMP>(PPC_SETTINGS_gasenin_l_djstra_omp);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = GaseninLDjstraOmpPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(DjkstraOmpPerf, GaseninLDjstraOmpPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace gasenin_l_djstra
