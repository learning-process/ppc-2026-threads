#include <gtest/gtest.h>

#include "dolov_v_crs_mat_mult_seq/common/include/common.hpp"
#include "dolov_v_crs_mat_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace dolov_v_crs_mat_mult_seq {

class DolovVCrsMatMultRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 100;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(DolovVCrsMatMultRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DolovVCrsMatMultSeq>(PPC_SETTINGS_dolov_v_crs_mat_mult_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DolovVCrsMatMultRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DolovVCrsMatMultRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dolov_v_crs_mat_mult_seq
