#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"
#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/omp/include/ops_omp.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "shkryleva_s_shell_sort_simple_merge/stl/include/ops_stl.hpp"
#include "shkryleva_s_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

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
    std::sort(expected_data_.begin(), expected_data_.end());
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

const std::array<TestType, 1> kPerfParam = {TestType{InType{}, OutType{}}};

const auto kPerfTasksSeq = ppc::util::AddFuncTask<ShkrylevaSShellMergeSEQ, InType>(
    kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kPerfTasksOmp = ppc::util::AddFuncTask<ShkrylevaSShellMergeOMP, InType>(
    kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kPerfTasksStl = ppc::util::AddFuncTask<ShkrylevaSShellMergeSTL, InType>(
    kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kPerfTasksTbb = ppc::util::AddFuncTask<ShkrylevaSShellMergeTBB, InType>(
    kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);
const auto kPerfTasksAll = ppc::util::AddFuncTask<ShkrylevaSShellMergeALL, InType>(
    kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);

std::vector<ppc::util::Task<InType, OutType>> kAllPerfTasks;
kAllPerfTasks.reserve(kPerfTasksSeq.size() + kPerfTasksOmp.size() + kPerfTasksStl.size() + kPerfTasksTbb.size() +
                      kPerfTasksAll.size());
kAllPerfTasks.insert(kAllPerfTasks.end(), kPerfTasksSeq.begin(), kPerfTasksSeq.end());
kAllPerfTasks.insert(kAllPerfTasks.end(), kPerfTasksOmp.begin(), kPerfTasksOmp.end());
kAllPerfTasks.insert(kAllPerfTasks.end(), kPerfTasksStl.begin(), kPerfTasksStl.end());
kAllPerfTasks.insert(kAllPerfTasks.end(), kPerfTasksTbb.begin(), kPerfTasksTbb.end());
kAllPerfTasks.insert(kAllPerfTasks.end(), kPerfTasksAll.begin(), kPerfTasksAll.end());

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = ShkrylevaSShellMergePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShkrylevaSShellMergePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
