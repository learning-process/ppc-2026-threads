#include <gtest/gtest.h>

#include "zhurin_i_gauss_kernel_seq/all/include/ops_all.hpp"
#include "zhurin_i_gauss_kernel_seq/common/include/common.hpp"
#include "zhurin_i_gauss_kernel_seq/omp/include/ops_omp.hpp"
#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"
#include "zhurin_i_gauss_kernel_seq/stl/include/ops_stl.hpp"
#include "zhurin_i_gauss_kernel_seq/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace zhurin_i_test_task_threads {

class ZhurinIRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(ZhurinIRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZhurinITestTaskALL, ZhurinITestTaskOMP, ZhurinITestTaskSEQ,
                                ZhurinITestTaskSTL, ZhurinITestTaskTBB>(PPC_SETTINGS_zhurin_i_gauss_kernel_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZhurinIRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZhurinIRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zhurin_i_test_task_threads
