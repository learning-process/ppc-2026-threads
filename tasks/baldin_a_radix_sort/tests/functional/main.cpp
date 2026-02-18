#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
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
    input_data_ = 
    return;
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

namespace {

TEST_P(BaldinARadixSortFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {  std::make_tuple("Empty", {}),
                                              std::make_tuple("Single_Element", {42}),
                                              std::make_tuple("Already_Sorted", {1, 2, 3, 4, 5, 10, 20}),
                                              std::make_tuple("Reverse_Sorted", {50, 40, 30, 20, 10, 5, 1}),
                                              std::make_tuple("Duplicates", {5, 1, 5, 2, 1, 5, 5}),
                                              std::make_tuple("Negative_Only", {-10, -50, -1, -100, -2}),
                                              std::make_tuple("Mixed_Sign", {-10, 50, -1, 0, 100, -200, 5}),
                                              std::make_tuple("Max_Min_Int", {INT_MAX, INT_MIN, 0, -1, 1})
                                          };

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<BaldinARadixSortSEQ, InType>(kTestParam, PPC_SETTINGS_baldin_a_radix_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BaldinARadixSortFuncTests::PrintFuncTestName<BaldinARadixSortFuncTests>;

INSTANTIATE_TEST_SUITE_P(BaldinARadixSortTests, BaldinARadixSortFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace baldin_a_radix_sort
