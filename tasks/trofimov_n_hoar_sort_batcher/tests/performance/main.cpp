#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"
#include "trofimov_n_hoar_sort_batcher/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace trofimov_n_hoar_sort_batcher {

class TrofimovNHoarSortBatcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const size_t size = 300000;

    std::mt19937 gen(12345);
    std::uniform_int_distribution<int> dist(-1000000, 1000000);

    input_data_.resize(size);
    for (auto &v : input_data_) {
      v = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(TrofimovNHoarSortBatcherPerfTests, RunPerformanceTests) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TrofimovNHoarSortBatcherSEQ>(PPC_SETTINGS_trofimov_n_hoar_sort_batcher);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TrofimovNHoarSortBatcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(TrofimovNHoarSortBatcherPerfTests, TrofimovNHoarSortBatcherPerfTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace trofimov_n_hoar_sort_batcher
