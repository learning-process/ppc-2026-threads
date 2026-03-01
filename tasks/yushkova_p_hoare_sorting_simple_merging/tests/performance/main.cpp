#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/seq/include/ops_seq.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

class YushkovaPRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    constexpr size_t kCount = 100000;
    input_data_.resize(kCount);

    std::mt19937 gen(2026);
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

    for (size_t i = 0; i < kCount; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == input_data_.size() && std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(YushkovaPRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, YushkovaPHoareSortingSimpleMergingSEQ>(
    PPC_SETTINGS_yushkova_p_hoare_sorting_simple_merging);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = YushkovaPRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, YushkovaPRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace yushkova_p_hoare_sorting_simple_merging
