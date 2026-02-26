#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1000000;  // 1 миллион элементов для теста производительности
  InType input_data_;

  void SetUp() override {
    input_data_.resize(kCount_);
    std::mt19937 gen(42);
    for (int i = 0; i < kCount_; ++i) {
      input_data_[i] = gen() % 100000;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(RadixSortPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RadixSortSimpleMergingSEQ>(
    "spichek_d_radix_sort_for_integers_with_simple_merging");

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(RunModeTests, RadixSortPerfTest, kGtestValues, RadixSortPerfTest::CustomPerfTestName);

}  // namespace

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
