#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>

#include "shekhirev_v_hoare_batcher_sort_seq/common/include/common.hpp"
#include "shekhirev_v_hoare_batcher_sort_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

class ShekhirevVRunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  ShekhirevVRunPerfTest() = default;

 protected:
  const size_t k_array_size = 200000;
  InType input_data;

  void SetUp() override {
    input_data.resize(k_array_size);
    for (size_t i = 0; i < k_array_size; ++i) {
      input_data[i] = static_cast<int>(k_array_size - i);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(ShekhirevVRunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ShekhirevHoareBatcherSortSEQ>(PPC_SETTINGS_shekhirev_v_hoare_batcher_sort_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShekhirevVRunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShekhirevVRunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shekhirev_v_hoare_batcher_sort_seq
