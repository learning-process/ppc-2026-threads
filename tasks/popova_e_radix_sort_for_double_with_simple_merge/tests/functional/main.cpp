#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "popova_e_radix_sort_for_double_with_simple_merge/common/include/common.hpp"
#include "popova_e_radix_sort_for_double_with_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace popova_e_radix_sort_for_double_with_simple_merge_threads {

class PopovaERunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == 1;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(PopovaERunFuncTestsSEQ, RadixSortSEQ) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    std::make_tuple(2, "2"),     std::make_tuple(3, "3"),       std::make_tuple(4, "4"),
    std::make_tuple(5, "5"),     std::make_tuple(10, "10_1"),   std::make_tuple(10, "10_2"),
    std::make_tuple(10, "10_3"), std::make_tuple(17, "17"),     std::make_tuple(50, "50"),
    std::make_tuple(100, "100"), std::make_tuple(1000, "1000"), std::make_tuple(5000, "5000")};

const auto kTestTasksList = ppc::util::AddFuncTask<PopovaERadixSorForDoubleWithSimpleMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_popova_e_radix_sort_for_double_with_simple_merge);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = PopovaERunFuncTestsSEQ::PrintFuncTestName<PopovaERunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(RadixSortSEQTests, PopovaERunFuncTestsSEQ, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace popova_e_radix_sort_for_double_with_simple_merge_threads
