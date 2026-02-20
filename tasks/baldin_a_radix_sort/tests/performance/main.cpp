#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "baldin_a_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace baldin_a_radix_sort {

class BaldinARadixSortPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int sz = 10'000'000;

    input_data_.resize(sz);

    std::mt19937 gen(42);
    std::uniform_int_distribution<int> dist(INT_MIN, INT_MAX);

    for (auto &val : input_data_) {
      val = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
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
