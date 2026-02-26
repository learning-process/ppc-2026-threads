#include <gtest/gtest.h>

// #include <algorithm>
#include <array>
#include <cstddef>
#include <string>
// #include <tuple>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergeFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = std::get<0>(param);
    expected_data_ = std::get<1>(param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.size(); i++) {
      if (output_data[i] != expected_data_[i]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

TEST_P(ShkrylevaSShellMergeFuncTests, radixTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    TestType{InType{1, 2, 3, 4, 5, 6}, OutType{1, 2, 3, 4, 5, 6}, std::string("Sorted")},
    TestType{InType{9, 8, 7, 6, 5, 4, 3, 2, 1}, OutType{1, 2, 3, 4, 5, 6, 7, 8, 9}, std::string("Reversed")},
    TestType{InType{1}, OutType{1}, std::string("One")},
    TestType{InType{2, 2, 2, 2, 2, 2}, OutType{2, 2, 2, 2, 2, 2}, std::string("Equal")},
    TestType{InType{2, 2, 44, 2, 3, 5, 1}, OutType{1, 2, 2, 2, 3, 5, 44}, std::string("EqualPlus")},
    TestType{InType{2, 1}, OutType{1, 2}, std::string("TwoElems")},
    TestType{InType{1, -2, 3, -5}, OutType{-5, -2, 1, 3}, std::string("Negative")},
    TestType{InType{1, 22, 13, 51, 2, 1, 2, 2, 34, 41}, OutType{1, 1, 2, 2, 2, 13, 22, 34, 41, 51},
             std::string("Mixed")}};

const auto kTestTasksList = ppc::util::AddFuncTask<ShkrylevaSShellMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_shkryleva_s_shell_sort_simple_merge);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = ShkrylevaSShellMergeFuncTests::PrintFuncTestName<ShkrylevaSShellMergeFuncTests>;

INSTANTIATE_TEST_SUITE_P(radixTest1, ShkrylevaSShellMergeFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
