#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>

#include "olesnitskiy_v_hoare_sort_simple_merge_seq/common/include/common.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge_seq {

namespace {
int NextPseudoRandom(int &state) {
  state = (state * 1103515245 + 12345) & 0x7fffffff;
  return state;
}
}  // namespace

class OlesnitskiyVRunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    constexpr size_t kSize = 100000;
    input_data_.resize(kSize);

    int state = 2026;
    std::ranges::generate(input_data_, [&]() {
      int value = NextPseudoRandom(state);
      return (value % 2000001) - 1000000;
    });
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(OlesnitskiyVRunPerfTestSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, OlesnitskiyVHoareSortSimpleMergeSEQ>(
    PPC_SETTINGS_olesnitskiy_v_hoare_sort_simple_merge_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = OlesnitskiyVRunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, OlesnitskiyVRunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace olesnitskiy_v_hoare_sort_simple_merge_seq
