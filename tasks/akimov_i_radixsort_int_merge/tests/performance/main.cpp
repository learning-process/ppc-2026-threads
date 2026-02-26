#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "akimov_i_radixsort_int_merge/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace akimov_i_radixsort_int_merge {

class AkimovIRadixSortIntMergePerfTests
    : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_;
  InType expected_sorted_;

  void SetUp() override {
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(-1000, 1000);
    input_data_.resize(kCount_);
    for (int &val : input_data_) {
      val = dist(gen);
    }
    expected_sorted_ = input_data_;
    std::sort(expected_sorted_.begin(), expected_sorted_.end());
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_sorted_;
  }

  InType GetTestInputData() final { return input_data_; }
};

TEST_P(AkimovIRadixSortIntMergePerfTests, RunPerfModes) { ExecuteTest(GetParam()); }

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, AkimovIRadixSortIntMergeSEQ>(
        PPC_SETTINGS_akimov_i_radixsort_int_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = AkimovIRadixSortIntMergePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, AkimovIRadixSortIntMergePerfTests,
                         kGtestValues, kPerfTestName);

}  // namespace

}  // namespace akimov_i_radixsort_int_merge