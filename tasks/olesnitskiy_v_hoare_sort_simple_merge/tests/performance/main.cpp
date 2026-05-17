#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge/omp/include/ops_omp.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge/seq/include/ops_seq.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

class OlesnitskiyVRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    constexpr std::size_t kCount = 100000;
    input_data_.resize(kCount);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

    for (std::size_t i = 0; i < kCount; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == input_data_.size() && std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(OlesnitskiyVRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, OlesnitskiyVHoareSortSimpleMergeSEQ, OlesnitskiyVHoareSortSimpleMergeOMP,
                                OlesnitskiyVHoareSortSimpleMergeTBB>(
        PPC_SETTINGS_olesnitskiy_v_hoare_sort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = OlesnitskiyVRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, OlesnitskiyVRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
