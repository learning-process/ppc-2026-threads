#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <ostream>
#include <string>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace ppc::util {
template <typename InType, typename OutType, typename TestType>
inline void PrintTo(const FuncTestParam<InType, OutType, TestType> &param, ::std::ostream *os) {
  *os << "FuncTestParam{"
      << "name=" << std::get<static_cast<std::size_t>(GTestParamIndex::kNameTest)>(param) << "}";
}
}  // namespace ppc::util

namespace sakharov_a_shell_sorting_with_merging_butcher {

class SakharovAShellButcherFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TestType MakeCase(const InType &input, const std::string &name) {
  OutType expected = input;
  std::ranges::sort(expected);
  return TestType{input, expected, name};
}

TEST_P(SakharovAShellButcherFuncTests, SortingWithButcherMerge) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    MakeCase({5, 1, 8, 3, 2, 7, 4, 6}, "basic_even"), MakeCase({9, -2, 4, -7, 0, 3, -1}, "negative_odd"),
    MakeCase({1, 1, 1, 1, 1, 1}, "duplicates"), MakeCase({42}, "single"), MakeCase({}, "empty")};

const auto kTestTasksList = ppc::util::AddFuncTask<SakharovAShellButcherSEQ, InType>(
    kTestParam, PPC_SETTINGS_sakharov_a_shell_sorting_with_merging_butcher);
const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kFuncTestName = SakharovAShellButcherFuncTests::PrintFuncTestName<SakharovAShellButcherFuncTests>;

INSTANTIATE_TEST_SUITE_P(ShellButcherSeqFunc, SakharovAShellButcherFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
