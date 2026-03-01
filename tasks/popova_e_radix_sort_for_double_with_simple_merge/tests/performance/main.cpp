#include <gtest/gtest.h>

// #include "popova_e_radix_sort_for_double_with_simple_merge/all/include/ops_all.hpp"
#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
// #include "popova_e_radix_sort_for_double_with_simple_merge/omp/include/ops_omp.hpp"
#include "popova_e_radix_sort_for_double_with_simple_merge/seq/include/ops_seq.hpp"
// #include "popova_e_radix_sort_for_double_with_simple_merge/stl/include/ops_stl.hpp"
// #include "popova_e_radix_sort_for_double_with_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

class PopovaERadixSortRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int kCount_ = 10000;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == 1;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PopovaERadixSortRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, PopovaERadixSorForDoubleWithSimpleMergeSEQ>(
    PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PopovaERadixSortRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PopovaERadixSortRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads