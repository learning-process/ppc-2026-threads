#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"
#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/omp/include/ops_omp.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "shkryleva_s_shell_sort_simple_merge/stl/include/ops_stl.hpp"
#include "shkryleva_s_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

// NOLINTNEXTLINE(misc-include-cleaner) — needed for std::vector via InType
#include <vector>

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    input_data_.resize(kCount_);
    expected_data_.resize(kCount_);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-100, 100);

    for (int i = 0; i < kCount_; ++i) {
      int number = dist(gen);
      input_data_[i] = number;
      expected_data_[i] = number;
    }
    std::ranges::sort(expected_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    return output_data == expected_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  const int kCount_ = 1000000;
  InType input_data_;
  OutType expected_data_;
};

TEST_P(ShkrylevaSShellMergePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

// --- Последовательная версия ---
const auto kPerfTasksSeq =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeSEQ>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
INSTANTIATE_TEST_SUITE_P(RunPerfSeq, ShkrylevaSShellMergePerfTests, ppc::util::TupleToGTestValues(kPerfTasksSeq),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- OpenMP версия ---
const auto kPerfTasksOmp =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeOMP>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
INSTANTIATE_TEST_SUITE_P(RunPerfOmp, ShkrylevaSShellMergePerfTests, ppc::util::TupleToGTestValues(kPerfTasksOmp),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- STL версия ---
const auto kPerfTasksStl =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeSTL>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
INSTANTIATE_TEST_SUITE_P(RunPerfStl, ShkrylevaSShellMergePerfTests, ppc::util::TupleToGTestValues(kPerfTasksStl),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- TBB версия ---
const auto kPerfTasksTbb =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeTBB>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
INSTANTIATE_TEST_SUITE_P(RunPerfTbb, ShkrylevaSShellMergePerfTests, ppc::util::TupleToGTestValues(kPerfTasksTbb),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- ALL версия ---
const auto kPerfTasksAll =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeALL>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
INSTANTIATE_TEST_SUITE_P(RunPerfAll, ShkrylevaSShellMergePerfTests, ppc::util::TupleToGTestValues(kPerfTasksAll),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

}  // namespace
}  // namespace shkryleva_s_shell_sort_simple_merge
