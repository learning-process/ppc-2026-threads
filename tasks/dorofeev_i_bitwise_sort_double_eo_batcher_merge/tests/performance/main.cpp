#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

class DorofeevIPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  // Вот эти две переменные, которые "потерял" редактор:
  const int k_count = 100000;
  InType input_data;

  void SetUp() override {
    // Честный случайный сид по стандарту безопасности
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-5000.0, 5000.0);

    input_data.resize(k_count);
    for (auto &val : input_data) {
      val = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // В перф тестах достаточно проверить, что массив просто отсортирован
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

namespace {

TEST_P(DorofeevIPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType,
                                                       /*DorofeevIBitwiseSortDoubleEOBatcherMergeALL,*/
                                                       /*DorofeevIBitwiseSortDoubleEOBatcherMergeOMP,*/
                                                       DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ
                                                       /*DorofeevIBitwiseSortDoubleEOBatcherMergeSTL,*/
                                                       /*DorofeevIBitwiseSortDoubleEOBatcherMergeTBB*/>(
    PPC_SETTINGS_dorofeev_i_bitwise_sort_double_eo_batcher_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = DorofeevIPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, DorofeevIPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
