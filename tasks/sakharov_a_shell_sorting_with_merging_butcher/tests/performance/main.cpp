#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ostream>

#include "sakharov_a_shell_sorting_with_merging_butcher/all/include/ops_all.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/omp/include/ops_omp.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/stl/include/ops_stl.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace ppc::util {
namespace {
template <typename InType, typename OutType>
void PrintTo(const PerfTestParam<InType, OutType> &param, ::std::ostream *os) {
  *os << "PerfTestParam{"
      << "name=" << std::get<static_cast<std::size_t>(GTestParamIndex::kNameTest)>(param) << "}";
}
}  // namespace
}  // namespace ppc::util

namespace sakharov_a_shell_sorting_with_merging_butcher {

class SakharovAShellButcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_output_;

  void SetUp() override {
    constexpr std::size_t kSize = 400000;
    std::uint32_t state = 12345;

    input_data_.resize(kSize);
    for (auto &value : input_data_) {
      state = (1664525U * state) + 1013904223U;
      value = static_cast<int>(state % 200001U) - 100000;
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

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SakharovAShellButcherSEQ, SakharovAShellButcherALL, SakharovAShellButcherOMP,
                                SakharovAShellButcherTBB, SakharovAShellButcherSTL>(
        PPC_SETTINGS_sakharov_a_shell_sorting_with_merging_butcher);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = SakharovAShellButcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ShellButcherSeqPerf, SakharovAShellButcherPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace sakharov_a_shell_sorting_with_merging_butcher
