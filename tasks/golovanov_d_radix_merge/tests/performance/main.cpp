#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>

// #include "golovanov_d_radix_merge/all/include/ops_all.hpp"
#include "golovanov_d_radix_merge/common/include/common.hpp"
// #include "golovanov_d_radix_merge/omp/include/ops_omp.hpp"
#include "golovanov_d_radix_merge/seq/include/ops_seq.hpp"
// #include "golovanov_d_radix_merge/stl/include/ops_stl.hpp"
// #include "golovanov_d_radix_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace golovanov_d_radix_merge {

class GolovanovDRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};
  int count_ = 500000;

  void SetUp() override {
    std::vector<double> data;
    data.reserve(2 * count_ + 1);
    for (int i = count_; i >= -count_; --i) {
      data.push_back(static_cast<double>(i));
    }

    input_data_ = std::move(data);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
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
