#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>

#include "ovchinnikov_m_shell_sort_batcher_merge_seq/common/include/common.hpp"
#include "ovchinnikov_m_shell_sort_batcher_merge_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge_seq {

class OvchinnikovMRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr std::size_t kSize = 100000;
    input_data_.resize(kSize);
    for (std::size_t i = 0; i < kSize; ++i) {
      input_data_[i] = static_cast<int>(kSize - i);
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == input_data_.size() && std::ranges::is_sorted(output_data.begin(), output_data.end());
  }

 private:
  InType input_data_;
};

TEST_P(OvchinnikovMRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, OvchinnikovMShellSortBatcherMergeSEQ>(
    PPC_SETTINGS_ovchinnikov_m_shell_sort_batcher_merge_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = OvchinnikovMRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ShellSortBatcherMergePerfTests, OvchinnikovMRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ovchinnikov_m_shell_sort_batcher_merge_seq
