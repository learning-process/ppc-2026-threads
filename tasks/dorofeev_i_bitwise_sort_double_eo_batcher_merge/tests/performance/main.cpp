#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

class DorofeevIPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected: 
  const int kCount_ = 100000; // Большой массив для проверки производительности
  InType input_data_{};

  void SetUp() override {
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-5000.0, 5000.0);
    
    input_data_.resize(kCount_);
    for (auto& val : input_data_) {
      val = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // В перф тестах достаточно проверить, что массив просто отсортирован
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

namespace {

TEST_P(DorofeevIPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, 
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