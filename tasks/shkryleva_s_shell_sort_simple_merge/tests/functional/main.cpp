#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"
#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/omp/include/ops_omp.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "shkryleva_s_shell_sort_simple_merge/stl/include/ops_stl.hpp"
#include "shkryleva_s_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

struct StaticTestCase {
  std::vector<int> input;
  std::string name;
};

const std::array<StaticTestCase, 8> kStaticTestCases = {
    StaticTestCase{.input = {1, 2, 3, 4, 5, 6}, .name = "sorted_asc"},
    StaticTestCase{.input = {9, 8, 7, 6, 5, 4, 3, 2, 1}, .name = "sorted_desc"},
    StaticTestCase{.input = {1}, .name = "single"},
    StaticTestCase{.input = {2, 2, 2, 2, 2, 2, 2, 2}, .name = "all_equal"},
    StaticTestCase{.input = {2, 2, 44, 2, 3, 5, 1}, .name = "random_with_duplicates"},
    StaticTestCase{.input = {2, 1}, .name = "two_unsorted"},
    StaticTestCase{.input = {1, -2, 3, -5}, .name = "mixed_signs"},
    StaticTestCase{.input = {1, 22, 13, 51, 2, 1, 2, 2, 34, 41}, .name = "longer_mixed"},
};

class ShkrylevaSShellMergeFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    size_t idx = std::get<0>(test_param);
    return kStaticTestCases.at(idx).name;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    size_t idx = std::get<0>(params);
    input_data_ = kStaticTestCases.at(idx).input;
    expected_data_ = input_data_;
    std::ranges::sort(expected_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

TEST_P(ShkrylevaSShellMergeFuncTests, shellMergeTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = []() {
  std::array<TestType, 8> arr;
  for (size_t i = 0; i < arr.size(); ++i) {
    arr[i] = std::make_tuple(i, kStaticTestCases[i].name);
  }
  return arr;
}();

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ShkrylevaSShellMergeSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge),
                                           ppc::util::AddFuncTask<ShkrylevaSShellMergeOMP, InType>(
                                               kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge),
                                           ppc::util::AddFuncTask<ShkrylevaSShellMergeSTL, InType>(
                                               kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge),
                                           ppc::util::AddFuncTask<ShkrylevaSShellMergeTBB, InType>(
                                               kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge),
                                           ppc::util::AddFuncTask<ShkrylevaSShellMergeALL, InType>(
                                               kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = ShkrylevaSShellMergeFuncTests::PrintFuncTestName<ShkrylevaSShellMergeFuncTests>;

INSTANTIATE_TEST_SUITE_P(shellMergeTests, ShkrylevaSShellMergeFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
