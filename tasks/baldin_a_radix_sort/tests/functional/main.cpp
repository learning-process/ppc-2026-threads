#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "baldin_a_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace baldin_a_radix_sort {

class BaldinARadixSortFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<1>(params);

    ref_data_ = input_data_;
    std::ranges::sort(ref_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == ref_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  InType ref_data_;
};

namespace {

TEST_P(BaldinARadixSortFuncTests, RadixSortTest) {
  ExecuteTest(GetParam());
}

std::vector<int> GenRandomVector(size_t size, int min_val, int max_val) {
  std::vector<int> res(size);

  const unsigned int a = 1103515245;
  const unsigned int c = 12345;
  unsigned int cur_val = 42;

  for (auto &val : res) {
    cur_val = (a * cur_val) + c;
    int range = max_val - min_val + 1;
    val = min_val + static_cast<int>(cur_val % static_cast<unsigned int>(range));
  }

  return res;
}

const std::array<TestType, 13> kTestParam = {
    std::make_tuple("Empty", std::vector<int>{}),
    std::make_tuple("Single_Element", std::vector<int>{42}),
    std::make_tuple("Already_Sorted", std::vector<int>{1, 2, 3, 4, 5, 10, 20}),
    std::make_tuple("Reverse_Sorted", std::vector<int>{50, 40, 30, 20, 10, 5, 1}),
    std::make_tuple("Duplicates", std::vector<int>{5, 1, 5, 2, 1, 5, 5}),
    std::make_tuple("Negative_Only", std::vector<int>{-10, -50, -1, -100, -2}),
    std::make_tuple("Mixed_Sign", std::vector<int>{-10, 50, -1, 0, 100, -200, 5}),
    std::make_tuple("Max_Min_Int", std::vector<int>{INT_MAX, INT_MIN, 0, -1, 1}),
    std::make_tuple("All_Zeros", std::vector<int>{0, 0, 0, 0, 0}),
    std::make_tuple("Alternating_Signs", std::vector<int>{1, -1, 2, -2, 3, -3}),
    std::make_tuple("Random_Small_Range", GenRandomVector(100, -10, 10)),
    std::make_tuple("Random_Large_Numbers", std::vector<int>{1000000, -1000000, 500000, -999999, 12345678}),
    std::make_tuple("Random_100_Elements", GenRandomVector(100, -1000, 1000))};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<BaldinARadixSortSEQ, InType>(kTestParam, PPC_SETTINGS_baldin_a_radix_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BaldinARadixSortFuncTests::PrintFuncTestName<BaldinARadixSortFuncTests>;

INSTANTIATE_TEST_SUITE_P(BaldinARadixSortTests, BaldinARadixSortFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace baldin_a_radix_sort
