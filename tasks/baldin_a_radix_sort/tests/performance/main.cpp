#include <gtest/gtest.h>

#include <algorithm>
#include <climits>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "baldin_a_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace baldin_a_radix_sort {

class BaldinARadixSortPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int sz = 10'000'000;

    input_data_.resize(sz);

    unsigned int cur_val = 42;

    for (auto &val : input_data_) {
      cur_val = (1664525 * cur_val) + 1013904223;
      val = static_cast<int>(cur_val);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(BaldinARadixSortPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, BaldinARadixSortSEQ>(PPC_SETTINGS_baldin_a_radix_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BaldinARadixSortPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BaldinARadixSortPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace baldin_a_radix_sort
