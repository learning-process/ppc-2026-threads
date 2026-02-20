#include <gtest/gtest.h>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/all/include/ops_all.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/omp/include/ops_omp.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/seq/include/ops_seq.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/stl/include/ops_stl.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

class DorofeevIRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(DorofeevIRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DorofeevIBitwiseSortDoubleEOBatcherMergeALL, DorofeevIBitwiseSortDoubleEOBatcherMergeOMP, DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ,
                                DorofeevIBitwiseSortDoubleEOBatcherMergeSTL, DorofeevIBitwiseSortDoubleEOBatcherMergeTBB>("dorofeev_i_bitwise_sort_double_eo_batcher_merge");

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DorofeevIRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DorofeevIRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
