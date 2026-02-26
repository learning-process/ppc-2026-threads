#include <gtest/gtest.h>

#include "golovanov_d_radix_merge/all/include/ops_all.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"
#include "golovanov_d_radix_merge/omp/include/ops_omp.hpp"
#include "golovanov_d_radix_merge/seq/include/ops_seq.hpp"
#include "golovanov_d_radix_merge/stl/include/ops_stl.hpp"
#include "golovanov_d_radix_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace golovanov_d_radix_merge {

class GolovanovDRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
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

TEST_P(GolovanovDRunPerfTestsThreads, RadixMergePerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, GolovanovDRadixMergeSEQ>(PPC_SETTINGS_golovanov_d_radix_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GolovanovDRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RadixMergePerf, GolovanovDRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace golovanov_d_radix_merge
