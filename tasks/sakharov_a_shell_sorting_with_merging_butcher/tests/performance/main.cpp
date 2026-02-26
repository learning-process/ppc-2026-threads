#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <random>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace ppc::util {
template <typename InType, typename OutType>
inline void PrintTo(const PerfTestParam<InType, OutType> &param, ::std::ostream *os) {
  *os << "PerfTestParam{"
      << "name=" << std::get<static_cast<std::size_t>(GTestParamIndex::kNameTest)>(param) << "}";
}
}  // namespace ppc::util

namespace sakharov_a_shell_sorting_with_merging_butcher {

class SakharovAShellButcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_output_;

  void SetUp() override {
    constexpr std::size_t kSize = 40000;
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<int> distribution(-100000, 100000);

    input_data_.resize(kSize);
    for (auto &value : input_data_) {
      value = distribution(generator);
    }

    expected_output_ = input_data_;
    std::ranges::sort(expected_output_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SakharovAShellButcherPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SakharovAShellButcherSEQ>(
    PPC_SETTINGS_sakharov_a_shell_sorting_with_merging_butcher);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = SakharovAShellButcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ShellButcherSeqPerf, SakharovAShellButcherPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace sakharov_a_shell_sorting_with_merging_butcher
