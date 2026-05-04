#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
// #include <ranges>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    input_data_.resize(kCount_);
    expected_data_.resize(kCount_);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-100, 100);

    for (int i = 0; i < kCount_; i++) {
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

    for (size_t i = 0; i < output_data.size(); i++) {
      if (output_data[i] != expected_data_[i]) {
        return false;
      }
    }

    return true;
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

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShkrylevaSShellMergeSEQ>(PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShkrylevaSShellMergePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShkrylevaSShellMergePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
