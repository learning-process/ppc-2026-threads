#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"
#include "papulina_y_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace papulina_y_radix_sort {

class PapulinaYRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 10000000;
  InType input_data_;
  std::vector<double> expected_result_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<double> dist(-1000.0, 1000.0);
    input_data_ = std::vector<double>(kCount_);
    for (int i = 0; i < kCount_; i++) {
      input_data_[i] = dist(gen);
    }
    expected_result_ = input_data_;
    std::ranges::sort(expected_result_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data == expected_result_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(PapulinaYRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PapulinaYRadixSortSEQ>(PPC_SETTINGS_papulina_y_radix_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PapulinaYRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PapulinaYRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace papulina_y_radix_sort
