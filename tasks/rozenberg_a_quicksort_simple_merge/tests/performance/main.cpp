#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "rozenberg_a_quicksort_simple_merge/all/include/ops_all.hpp"
#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"
#include "rozenberg_a_quicksort_simple_merge/omp/include/ops_omp.hpp"
#include "rozenberg_a_quicksort_simple_merge/seq/include/ops_seq.hpp"
#include "rozenberg_a_quicksort_simple_merge/stl/include/ops_stl.hpp"
#include "rozenberg_a_quicksort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace rozenberg_a_quicksort_simple_merge {

class RozenbergARunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    input_data_.clear();
    output_data_.clear();

    constexpr int kSize = 10000000;
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(-10000000, 10000000);

    InType input_data(kSize);
    for (int i = 0; i < kSize; i++) {
      input_data[i] = dist(rng);
    }

    OutType output_data = input_data;
    std::ranges::sort(output_data);

    input_data_ = input_data;
    output_data_ = output_data;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (ppc::util::IsUnderMpirun() && ppc::util::GetMPIRank() != 0) {
      return true;
    }
    return (output_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType output_data_;
};

TEST_P(RozenbergARunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RozenbergAQuicksortSimpleMergeSEQ, RozenbergAQuicksortSimpleMergeOMP,
                                RozenbergAQuicksortSimpleMergeTBB, RozenbergAQuicksortSimpleMergeSTL,
                                RozenbergAQuicksortSimpleMergeALL>(PPC_SETTINGS_rozenberg_a_quicksort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RozenbergARunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RozenbergARunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace rozenberg_a_quicksort_simple_merge
