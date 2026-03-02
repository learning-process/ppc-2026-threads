#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "likhanov_m_hoare_batcher_sort/common/include/common.hpp"
#include "likhanov_m_hoare_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace likhanov_m_hoare_batcher_sort {

class LikhanovMHoareBatcherSortRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_count = 2000;
  const int k_size = 1000000;

  InType input_data;
  OutType expected_data;

  void SetUp() override {
    input_data.resize(k_size);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

    for (auto &v : input_data) {
      v = dist(gen);
    }

    expected_data = input_data;
    std::ranges::sort(expected_data);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_data;
  }

  InType GetTestInputData() final {
    return input_data;
  }

  [[nodiscard]] int GetRepeatCount() const {
    return k_count;
  }
};

TEST_P(LikhanovMHoareBatcherSortRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LikhanovMHoareBatcherSortSEQ>(PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LikhanovMHoareBatcherSortRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, LikhanovMHoareBatcherSortRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace likhanov_m_hoare_batcher_sort
