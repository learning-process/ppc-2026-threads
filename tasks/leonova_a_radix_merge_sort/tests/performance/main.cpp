#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "leonova_a_radix_merge_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace leonova_a_radix_merge_sort {

namespace {  // namespace

std::vector<int64_t> GenerateVector(size_t size) {
  std::vector<int64_t> result(size);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(-1000000, 1000000);

  for (size_t i = 0; i < size; ++i) {
    result[i] = dist(gen);
  }
  return result;
}

}  // namespace

class LeonovaARunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  std::vector<TestType> test_data_;
  size_t current_index_ = 0;

  void SetUp() override {
    std::vector<size_t> sizes = {10000, 100000, 1000000, 10000000};

    for (size_t size : sizes) {
      auto input = GenerateVector(size);
      auto expected = input;
      std::ranges::sort(expected);
      test_data_.emplace_back(input, expected);
    }

    current_index_ = 0;
  }

  bool CheckTestOutputData(OutType &output) final {
    if (current_index_ == 0) {
      return false;
    }

    const auto &expected = std::get<1>(test_data_[current_index_ - 1]);

    if (output.size() != expected.size()) {
      return false;
    }

    for (size_t i = 1; i < output.size(); ++i) {
      if (output[i - 1] > output[i]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    if (current_index_ < test_data_.size()) {
      return std::get<0>(test_data_[current_index_++]);
    }
    return std::get<0>(test_data_.back());
  }
};

TEST_P(LeonovaARunPerfTests, RunPerfRadix) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LeonovaARadixMergeSortSEQ>(PPC_SETTINGS_leonova_a_radix_merge_sort);

INSTANTIATE_TEST_SUITE_P(RunPerfRadixTests, LeonovaARunPerfTests, ppc::util::TupleToGTestValues(kAllPerfTasks),
                         LeonovaARunPerfTests::CustomPerfTestName);
}  // namespace
}  // namespace leonova_a_radix_merge_sort
