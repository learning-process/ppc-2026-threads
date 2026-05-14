#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "olesnitskiy_v_hoare_sort_simple_merge/common/include/common.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge/omp/include/ops_omp.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge {

class OlesnitskiyVRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected_data = input_data_;
    std::ranges::sort(expected_data);
    return output_data == expected_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(OlesnitskiyVRunFuncTestsThreads, HoareSortSimpleMerging) {
  ExecuteTest(GetParam());
}

TestType MakeTest(InType input, std::string name) {
  return {std::move(input), std::move(name)};
}

const std::array<TestType, 15> kTestParam = {
    MakeTest(std::vector<int>{42}, "single"),
    MakeTest(std::vector<int>{2, 1}, "two_elements"),
    MakeTest(std::vector<int>{1, 2, 3, 4, 5}, "already_sorted"),
    MakeTest(std::vector<int>{5, 4, 3, 2, 1}, "reverse_sorted"),
    MakeTest(std::vector<int>{5, 1, 5, 1, 5, 1}, "duplicates"),
    MakeTest(std::vector<int>{0, -1, 7, -5, 2, -3}, "mixed_signs"),
    MakeTest(std::vector<int>{10, 3, 8, 6, 4, 9, 2, 7, 1, 5}, "random_10"),
    MakeTest(std::vector<int>{100, 1, 50, 2, 75, 3, 60, 4, 20, 5, 30}, "odd_count"),
    MakeTest(std::vector<int>{9, 9, 8, 8, 7, 7, 6, 6, 5, 5}, "pair_duplicates"),
    MakeTest(std::vector<int>{1000, -1000, 500, -500, 0, 250, -250}, "wide_range"),
    MakeTest(std::vector<int>(64, 7), "one_block_equal"),
    MakeTest(std::vector<int>{64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43,
                              42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21,
                              20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1,  0},
             "crosses_block"),
    MakeTest(std::vector<int>{std::numeric_limits<int>::max(), 0, std::numeric_limits<int>::min(), -1, 1},
             "int_limits"),
    MakeTest(std::vector<int>{3, 3, 3, 2, 2, 1, 1, 0, 0, -1, -1}, "many_equal_runs"),
    MakeTest(std::vector<int>{17, -4, 23, 0, 17, -4, 99, -100, 8, 8, 42, -100, 5}, "repeated_mixed")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<OlesnitskiyVHoareSortSimpleMergeSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_olesnitskiy_v_hoare_sort_simple_merge),
                                           ppc::util::AddFuncTask<OlesnitskiyVHoareSortSimpleMergeOMP, InType>(
                                               kTestParam, PPC_SETTINGS_olesnitskiy_v_hoare_sort_simple_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = OlesnitskiyVRunFuncTestsThreads::PrintFuncTestName<OlesnitskiyVRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(HoareSortSimpleMergingTests, OlesnitskiyVRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace olesnitskiy_v_hoare_sort_simple_merge
