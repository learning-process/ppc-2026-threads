/*
#include <gtest/gtest.h>

// #include "barkalova_m_mult_matrix_ccs/all/include/ops_all.hpp"
#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
// #include "barkalova_m_mult_matrix_ccs/omp/include/ops_omp.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
// #include "barkalova_m_mult_matrix_ccs/stl/include/ops_stl.hpp"
// #include "barkalova_m_mult_matrix_ccs/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMMultMatrixCcsPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
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

TEST_P(BarkalovaMMultMatrixCcsPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BarkalobaMMultMatrixCcsSEQ>(PPC_SETTINGS_barkalova_m_mult_matrix_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BarkalovaMMultMatrixCcsPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BarkalovaMMultMatrixCcsPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace barkalova_m_mult_matrix_ccs
*/



