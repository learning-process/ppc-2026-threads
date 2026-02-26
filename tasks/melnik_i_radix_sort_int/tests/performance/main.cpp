#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "melnik_i_radix_sort_int/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace melnik_i_radix_sort_int {

namespace {

constexpr std::size_t kArraySize = 1'000'000;
constexpr unsigned int kSeed = 42;

std::vector<int> GenerateRandomData(std::size_t size, unsigned int seed) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(-1000000, 1000000);
  std::vector<int> data(size);
  for (auto &val : data) {
    val = dist(gen);
  }
  return data;
}

}  // namespace

class MelnikIRadixSortIntPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    input_data_ = GenerateRandomData(kArraySize, kSeed);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(MelnikIRadixSortIntPerfTests, SortPerformance) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MelnikIRadixSortIntSEQ>(PPC_SETTINGS_melnik_i_radix_sort_int);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MelnikIRadixSortIntPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RadixSortPerfTests, MelnikIRadixSortIntPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace melnik_i_radix_sort_int
