#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "sosnina_a_radix_simple_merge/common/include/common.hpp"
#include "sosnina_a_radix_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sosnina_a_radix_simple_merge {

struct StaticTestCase {
  std::vector<int> input;
  std::string name;
};

const std::array<StaticTestCase, 48> kStaticTestCases = {
    StaticTestCase{.input = std::vector<int>{1}, .name = "single"},
    StaticTestCase{.input = std::vector<int>{2, 1}, .name = "two_unsorted"},
    StaticTestCase{.input = std::vector<int>{1, 2}, .name = "two_sorted"},
    StaticTestCase{.input = std::vector<int>{5, 3, 8, 2, 1}, .name = "small_mixed"},
    StaticTestCase{.input = std::vector<int>{1, 2, 3, 4, 5}, .name = "sorted_5"},
    StaticTestCase{.input = std::vector<int>{5, 4, 3, 2, 1}, .name = "reversed_5"},
    StaticTestCase{.input = std::vector<int>{7, 7, 7, 7}, .name = "all_same"},
    StaticTestCase{.input = std::vector<int>{0, 0, 0}, .name = "all_zero"},
    StaticTestCase{.input = std::vector<int>{100, -50, 0, 25, -25}, .name = "with_zero"},
    StaticTestCase{.input = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, .name = "sorted_10"},
    StaticTestCase{.input = std::vector<int>{10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, .name = "reversed_10"},
    StaticTestCase{.input = std::vector<int>{3, 1, 4, 1, 5, 9, 2, 6}, .name = "duplicates"},
    StaticTestCase{.input = std::vector<int>{-1, -2, -3, -4, -5}, .name = "negative_sorted"},
    StaticTestCase{.input = std::vector<int>{-5, -4, -3, -2, -1}, .name = "negative_reversed"},
    StaticTestCase{.input = std::vector<int>{42, 42, 42, 42, 42, 42}, .name = "all_42"},
    StaticTestCase{.input = std::vector<int>{INT_MAX}, .name = "max_int"},
    StaticTestCase{.input = std::vector<int>{INT_MIN}, .name = "min_int"},
    StaticTestCase{.input = std::vector<int>{INT_MIN, INT_MAX}, .name = "min_max"},
    StaticTestCase{.input = std::vector<int>{INT_MAX, INT_MIN}, .name = "max_min"},
    StaticTestCase{.input = std::vector<int>{0, INT_MAX, INT_MIN, 1, -1}, .name = "extreme_mixed"},
    StaticTestCase{.input = std::vector<int>{1, 2, 2, 3, 3, 3, 4, 4, 4, 4}, .name = "ascending_dups"},
    StaticTestCase{.input = std::vector<int>{4, 4, 4, 4, 3, 3, 3, 2, 2, 1}, .name = "descending_dups"},
    StaticTestCase{.input = std::vector<int>{17, 23, 5, 11, 29, 7, 13, 19, 31, 2}, .name = "primes"},
    StaticTestCase{.input = std::vector<int>{256, 128, 64, 32, 16, 8, 4, 2, 1}, .name = "powers_of_2"},
    StaticTestCase{.input = std::vector<int>{1, 2, 4, 8, 16, 32, 64, 128, 256}, .name = "powers_sorted"},
    StaticTestCase{.input = std::vector<int>{-100, -50, 0, 50, 100}, .name = "symmetric"},
    StaticTestCase{.input = std::vector<int>{1, 1, 2, 2, 3, 3}, .name = "pairs"},
    StaticTestCase{.input = std::vector<int>{3, 2, 1, 1, 2, 3}, .name = "palindrome_unsorted"},
    StaticTestCase{.input = std::vector<int>{1, 1, 1, 2, 2, 2, 3, 3, 3}, .name = "triples"},
    StaticTestCase{.input = std::vector<int>{9, 7, 5, 3, 1, 2, 4, 6, 8}, .name = "odd_even"},
    StaticTestCase{.input = std::vector<int>{15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, .name = "reversed_15"},
    StaticTestCase{.input = std::vector<int>{1, 3, 5, 7, 9, 2, 4, 6, 8, 10}, .name = "interleaved"},
    StaticTestCase{.input = std::vector<int>{100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90}, .name = "decreasing_11"},
    StaticTestCase{.input = std::vector<int>{10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, .name = "increasing_10"},
    StaticTestCase{.input = std::vector<int>{5, 10, 15, 20, 25, 30, 35, 40, 45, 50}, .name = "step_5"},
    StaticTestCase{.input = std::vector<int>{-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, .name = "even_range"},
    StaticTestCase{.input = std::vector<int>{1, 1, 1, 1, 1, 2, 2, 2, 2, 2}, .name = "two_groups"},
    StaticTestCase{.input = std::vector<int>{2, 1, 2, 1, 2, 1, 2, 1}, .name = "alternating"},
    StaticTestCase{.input = std::vector<int>{1000, 500, 250, 125, 62, 31, 15, 7, 3, 1}, .name = "halving"},
    StaticTestCase{.input = std::vector<int>{1, 10, 100, 1000, 10000}, .name = "powers_of_10"},
    StaticTestCase{.input = std::vector<int>{10000, 1000, 100, 10, 1}, .name = "powers_of_10_rev"},
    StaticTestCase{.input = std::vector<int>{0, 1, 0, 1, 0, 1, 0, 1, 0, 1}, .name = "zero_one_alt"},
    StaticTestCase{.input = std::vector<int>{-5, 5, -4, 4, -3, 3, -2, 2, -1, 1}, .name = "pos_neg_alt"},
    StaticTestCase{.input = std::vector<int>{777, 777, 777}, .name = "triple_777"},
    StaticTestCase{.input = std::vector<int>{1, 2, 3, 2, 1}, .name = "mountain"},
    StaticTestCase{.input = std::vector<int>{3, 2, 1, 2, 3}, .name = "valley"},
    StaticTestCase{.input = std::vector<int>{6, 2, 9, 1, 5, 3, 8, 4, 7}, .name = "random_9"},
    StaticTestCase{.input = std::vector<int>{12, 5, 18, 3, 9, 14, 7, 21, 6, 11}, .name = "random_10"},
};

class SosninaARunFuncTestsRadixSort : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    size_t idx = std::get<0>(params);
    input_data_ = kStaticTestCases.at(idx).input;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    if (!std::ranges::is_sorted(output_data)) {
      return false;
    }
    std::vector<int> expected = input_data_;
    std::ranges::sort(expected);
    return output_data == expected;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(SosninaARunFuncTestsRadixSort, RadixSortSimpleMerge) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 48> kTestParam = []() {
  std::array<TestType, 48> arr;
  for (size_t i = 0; i < 48; ++i) {
    arr.at(i) = std::make_tuple(i, kStaticTestCases.at(i).name);
  }
  return arr;
}();

const auto kTestTasksList =
    ppc::util::AddFuncTask<SosninaATestTaskSEQ, InType>(kTestParam, PPC_SETTINGS_sosnina_a_radix_simple_merge);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SosninaARunFuncTestsRadixSort::PrintFuncTestName<SosninaARunFuncTestsRadixSort>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, SosninaARunFuncTestsRadixSort, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sosnina_a_radix_simple_merge
