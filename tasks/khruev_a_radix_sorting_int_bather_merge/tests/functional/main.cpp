#include <gtest/gtest.h>

// #include <algorithm>
#include <array>
#include <cstddef>
#include <string>
// #include <tuple>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"
#include "khruev_a_radix_sorting_int_bather_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

class KhruevARadixSortingIntBatherMergeFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(KhruevARadixSortingIntBatherMergeFuncTests, radixTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    TestType{InType{1, 2, 3, 4, 5, 6}, OutType{1, 2, 3, 4, 5, 6}, "Sorted"},
    TestType{InType{9, 8, 7, 6, 5, 4, 3, 2, 1}, OutType{1, 2, 3, 4, 5, 6, 7, 8, 9}, "Reversed"},
    TestType{InType{1}, OutType{1}, "One"},
    TestType{InType{2, 2, 2, 2, 2, 2}, OutType{2, 2, 2, 2, 2, 2}, "odinakovie"},
    TestType{InType{2, 2, 44, 2, 3, 5, 1}, OutType{1, 2, 2, 2, 3, 5, 44}, "odinakovie_plus"},
    TestType{InType{2, 1}, OutType{1, 2}, "two_elems"},
    TestType{InType{1, -2, 3, -5}, OutType{-5, -2, 1, 3}, "negative"},
    TestType{InType{1, 22, 13, 51, 2, 1, 2, 2, 34, 41}, OutType{1, 1, 2, 2, 2, 13, 22, 34, 41, 51}, "raznie"}

};

const auto kTestTasksList = ppc::util::AddFuncTask<KhruevARadixSortingIntBatherMergeSEQ, InType>(
    kTestParam, PPC_SETTINGS_khruev_a_radix_sorting_int_bather_merge);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName =
    KhruevARadixSortingIntBatherMergeFuncTests::PrintFuncTestName<KhruevARadixSortingIntBatherMergeFuncTests>;

INSTANTIATE_TEST_SUITE_P(radixTest1, KhruevARadixSortingIntBatherMergeFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace khruev_a_radix_sorting_int_bather_merge
