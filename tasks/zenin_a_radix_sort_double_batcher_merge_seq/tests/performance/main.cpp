#include <gtest/gtest.h>

//#include "zenin_a_radix_sort_double_batcher_merge_seq/all/include/ops_all.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_seq/common/include/common.hpp"
//#include "zenin_a_radix_sort_double_batcher_merge_seq/omp/include/ops_omp.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_seq/seq/include/ops_seq.hpp"
//#include "zenin_a_radix_sort_double_batcher_merge_seq/stl/include/ops_stl.hpp"
//#include "zenin_a_radix_sort_double_batcher_merge_seq/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

class ZeninARadixSortDoubleBatcherMergePerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(ZeninARadixSortDoubleBatcherMergePerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZeninARadixSortDoubleBatcherMerge_SEQSEQ>(PPC_SETTINGS_zenin_a_radix_sort_double_batcher_merge_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZeninARadixSortDoubleBatcherMergePerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZeninARadixSortDoubleBatcherMergePerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
