#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <random>

#include "shelenkova_m_shell_sort_simple_merge/common/include/common.hpp"
#include "shelenkova_m_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shelenkova_m_shell_sort_simple_merge {

namespace {

constexpr double kNanosToSeconds = 1e-9;
constexpr int kMinRandomValue = -1000000;
constexpr int kMaxRandomValue = 1000000;

}  // namespace

class ShelenkovaMRunPerfTestShellSort : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 100000;
  InType input_data_;

 protected:
  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attr) override {
    const auto t0 = std::chrono::steady_clock::now();
    perf_attr.current_timer = [t0] {
      auto now = std::chrono::steady_clock::now();
      auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - t0).count();
      return static_cast<double>(ns) * kNanosToSeconds;
    };
    perf_attr.num_running = 1;
  }

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(kMinRandomValue, kMaxRandomValue);

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

TEST_P(ShelenkovaMRunPerfTestShellSort, RunPerfShellSort) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, ShelenkovaMShellSortSimpleMergeSEQ>(
    PPC_SETTINGS_shelenkova_m_shell_sort_simple_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShelenkovaMRunPerfTestShellSort::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(ShellSortPerfTests, ShelenkovaMRunPerfTestShellSort, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shelenkova_m_shell_sort_simple_merge
