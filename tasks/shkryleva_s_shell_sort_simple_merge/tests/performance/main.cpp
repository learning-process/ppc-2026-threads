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
#include "util/include/func_test_util.hpp"
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

// Один параметр (данные будут сгенерированы в SetUp)
const std::array<TestType, 1> kPerfParam = {TestType{InType{}, OutType{}}};

// --- Последовательная версия ---
INSTANTIATE_TEST_SUITE_P(RunPerfSeq, ShkrylevaSShellMergePerfTests,
                         ppc::util::TupleToGTestValues(ppc::util::AddFuncTask<ShkrylevaSShellMergeSEQ, InType>(
                             kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge)),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- OpenMP версия ---
INSTANTIATE_TEST_SUITE_P(RunPerfOmp, ShkrylevaSShellMergePerfTests,
                         ppc::util::TupleToGTestValues(ppc::util::AddFuncTask<ShkrylevaSShellMergeOMP, InType>(
                             kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge)),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- STL версия ---
INSTANTIATE_TEST_SUITE_P(RunPerfStl, ShkrylevaSShellMergePerfTests,
                         ppc::util::TupleToGTestValues(ppc::util::AddFuncTask<ShkrylevaSShellMergeSTL, InType>(
                             kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge)),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- TBB версия ---
INSTANTIATE_TEST_SUITE_P(RunPerfTbb, ShkrylevaSShellMergePerfTests,
                         ppc::util::TupleToGTestValues(ppc::util::AddFuncTask<ShkrylevaSShellMergeTBB, InType>(
                             kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge)),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

// --- ALL версия ---
INSTANTIATE_TEST_SUITE_P(RunPerfAll, ShkrylevaSShellMergePerfTests,
                         ppc::util::TupleToGTestValues(ppc::util::AddFuncTask<ShkrylevaSShellMergeALL, InType>(
                             kPerfParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge)),
                         ShkrylevaSShellMergePerfTests::CustomPerfTestName);

}  // namespace
}  // namespace shkryleva_s_shell_sort_simple_merge
