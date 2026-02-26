#include <gtest/gtest.h>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "akimov_i_radixsort_int_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace akimov_i_radixsort_int_merge {

class ExampleRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(ExampleRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, AkimovIRadixSortIntMergeSEQ>(PPC_SETTINGS_akimov_i_radixsort_int_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ExampleRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ExampleRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace akimov_i_radixsort_int_merge
