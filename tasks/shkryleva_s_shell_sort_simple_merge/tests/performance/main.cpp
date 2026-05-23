#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"
#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/omp/include/ops_omp.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "shkryleva_s_shell_sort_simple_merge/stl/include/ops_stl.hpp"
#include "shkryleva_s_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    input_data_.resize(kCount_);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

    for (size_t i = 0; i < kCount_; ++i) {
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

 private:
  static constexpr size_t kCount_ = 100000;
  InType input_data_;
};

TEST_P(ShkrylevaSShellMergePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeALL, ShkrylevaSShellMergeOMP, ShkrylevaSShellMergeSEQ,
                                ShkrylevaSShellMergeSTL, ShkrylevaSShellMergeTBB>(
        PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = ShkrylevaSShellMergePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShkrylevaSShellMergePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
