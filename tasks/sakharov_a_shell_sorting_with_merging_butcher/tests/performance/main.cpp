#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

class ShellButcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_output_;

  void SetUp() override {
    constexpr size_t kSize = 40000;
    std::mt19937 generator(2026);
    std::uniform_int_distribution<int> distribution(-100000, 100000);

    input_data_.resize(kSize);
    for (auto &value : input_data_) {
      value = distribution(generator);
    }

    expected_output_ = input_data_;
    std::sort(expected_output_.begin(), expected_output_.end());
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ShellButcherPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SakharovAShellButcherSEQ>(
    PPC_SETTINGS_sakharov_a_shell_sorting_with_merging_butcher);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = ShellButcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ShellButcherSeqPerf, ShellButcherPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace sakharov_a_shell_sorting_with_merging_butcher
