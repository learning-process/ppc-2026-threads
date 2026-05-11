#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>

#include "olesnitskiy_v_hoare_sort_simple_merge_seq/common/include/common.hpp"
#include "olesnitskiy_v_hoare_sort_simple_merge_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace olesnitskiy_v_hoare_sort_simple_merge_seq {

class OlesnitskiyVRunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(param);
    expected_data_ = std::get<1>(param);
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

TEST_P(OlesnitskiyVRunFuncTestsSEQ, HoareSortSimpleMergeSeq) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    TestType{InType{5}, OutType{5}, "single_element"},
    TestType{InType{2, 1}, OutType{1, 2}, "two_elements"},
    TestType{InType{1, 2, 3, 4, 5}, OutType{1, 2, 3, 4, 5}, "already_sorted"},
    TestType{InType{9, 7, 5, 3, 1}, OutType{1, 3, 5, 7, 9}, "reverse_sorted"},
    TestType{InType{4, 1, 3, 4, 2, 3, 1}, OutType{1, 1, 2, 3, 3, 4, 4}, "with_duplicates"},
    TestType{InType{-3, 0, 5, -10, 8, -1}, OutType{-10, -3, -1, 0, 5, 8}, "negative_values"},
    TestType{InType{10, -5, 0, 10, -5, 8, 2, 2}, OutType{-5, -5, 0, 2, 2, 8, 10, 10}, "mixed_values"},
    TestType{InType{15, 3, 27, 9, 1, 8, 2, 11, 6, 4, 5}, OutType{1, 2, 3, 4, 5, 6, 8, 9, 11, 15, 27}, "odd_size"}};

const auto kTestTasksList = ppc::util::AddFuncTask<OlesnitskiyVHoareSortSimpleMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_olesnitskiy_v_hoare_sort_simple_merge_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = OlesnitskiyVRunFuncTestsSEQ::PrintFuncTestName<OlesnitskiyVRunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(HoareSimpleMergeFuncTests, OlesnitskiyVRunFuncTestsSEQ, kGtestValues, kTestName);

}  // namespace

}  // namespace olesnitskiy_v_hoare_sort_simple_merge_seq
