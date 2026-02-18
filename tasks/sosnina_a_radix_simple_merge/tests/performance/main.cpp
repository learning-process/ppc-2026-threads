#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

#include "sosnina_a_radix_simple_merge/common/include/common.hpp"
#include "sosnina_a_radix_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sosnina_a_radix_simple_merge {

class SosninaARunPerfTestRadixSort : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 100000;
  InType input_data_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-100000, 100000);

    input_data_.resize(kCount);
    for (size_t i = 0; i < kCount; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SosninaARunPerfTestRadixSort, RunPerfRadixSort) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SosninaATestTaskSEQ>(PPC_SETTINGS_sosnina_a_radix_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SosninaARunPerfTestRadixSort::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RadixSortPerfTests, SosninaARunPerfTestRadixSort, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sosnina_a_radix_simple_merge
