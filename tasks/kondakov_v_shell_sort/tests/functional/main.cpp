#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "kondakov_v_shell_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kondakov_v_shell_sort {

struct StaticTestCase {
  std::vector<int> input;
  std::string name;
};

const std::array<StaticTestCase, 12> kStaticTestCases = {
    StaticTestCase{.input = std::vector<int>{1}, .name = "single"},
    StaticTestCase{.input = std::vector<int>{2, 1}, .name = "two_unsorted"},
    StaticTestCase{.input = std::vector<int>{1, 2, 3, 4, 5}, .name = "already_sorted"},
    StaticTestCase{.input = std::vector<int>{5, 4, 3, 2, 1}, .name = "reverse_sorted"},
    StaticTestCase{.input = std::vector<int>{5, 3, 8, 2, 1}, .name = "small_random"},
    StaticTestCase{.input = std::vector<int>{7, 7, 7, 7, 7}, .name = "all_equal"},
    StaticTestCase{.input = std::vector<int>{-5, -1, -10, -3, -2}, .name = "all_negative"},
    StaticTestCase{.input = std::vector<int>{-3, 0, 5, -1, 2, -2}, .name = "mixed_signs"},
    StaticTestCase{.input = std::vector<int>{10, -10, 10, -10, 0}, .name = "duplicates_with_zero"},
    StaticTestCase{.input = std::vector<int>{9, 1, 8, 2, 7, 3, 6, 4, 5}, .name = "zigzag"},
    StaticTestCase{.input = std::vector<int>{1000, 500, 250, 125, 62, 31, 15, 7, 3, 1}, .name = "halving"},
    StaticTestCase{.input = std::vector<int>{42, -42, 17, 0, -1, 999, -999, 5}, .name = "extreme_mix"},
};

class KondakovVRunFuncTestsShellSort : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const size_t idx = std::get<0>(params);
    input_data_ = kStaticTestCases.at(idx).input;

    expected_result_ = input_data_;
    std::ranges::sort(expected_result_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_result_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_result_;
};

namespace {

TEST_P(KondakovVRunFuncTestsShellSort, ShellSortSimpleMerge) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = []() {
  std::array<TestType, 12> arr;
  for (size_t i = 0; i < arr.size(); ++i) {
    arr.at(i) = std::make_tuple(i, kStaticTestCases.at(i).name);
  }
  return arr;
}();

const auto kTestTasksList =
    ppc::util::AddFuncTask<KondakovVShellSortSEQ, InType>(kTestParam, PPC_SETTINGS_kondakov_v_shell_sort);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KondakovVRunFuncTestsShellSort::PrintFuncTestName<KondakovVRunFuncTestsShellSort>;

INSTANTIATE_TEST_SUITE_P(ShellSortTests, KondakovVRunFuncTestsShellSort, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kondakov_v_shell_sort
