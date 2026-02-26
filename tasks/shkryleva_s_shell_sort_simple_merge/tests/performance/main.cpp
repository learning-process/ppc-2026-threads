#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
#include <ranges>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_output_;

  void SetUp() override {
    constexpr size_t kSize = 40000;
    std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(-100000, 100000);

    input_data_.resize(kSize);
    for (auto &value : input_data_) {
      value = distribution(generator);
    }

    expected_output_ = input_data_;
    std::ranges::sort(expected_output_.begin(), expected_output_.end());
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ShkrylevaSShellMergePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeSEQ>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = ShkrylevaSShellMergePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ShellMergeSeqPerf, ShkrylevaSShellMergePerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace shkryleva_s_shell_sort_simple_merge
