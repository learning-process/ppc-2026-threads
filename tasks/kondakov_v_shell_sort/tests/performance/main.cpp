#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "kondakov_v_shell_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kondakov_v_shell_sort {

class KondakovVRunPerfTestsShellSort : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 5000000;
  InType input_data_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

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

TEST_P(KondakovVRunPerfTestsShellSort, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KondakovVShellSortSEQ>(PPC_SETTINGS_kondakov_v_shell_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KondakovVRunPerfTestsShellSort::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KondakovVRunPerfTestsShellSort, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kondakov_v_shell_sort
