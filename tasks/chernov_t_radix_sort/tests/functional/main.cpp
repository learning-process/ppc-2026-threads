#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "chernov_t_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(ChernovTRadixSortFuncTests, RadixSortTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParam = {
    std::make_tuple("NoElements", std::vector<int>{}),
    std::make_tuple("JustOneItem", std::vector<int>{42}),
    std::make_tuple("AscendingOrder", std::vector<int>{1, 2, 3, 4, 5, 10, 20}),
    std::make_tuple("OnlyNegatives", std::vector<int>{-10, -50, -1, -100, -2}),
    std::make_tuple("PosAndNegMixed", std::vector<int>{-10, 50, -1, 0, 100, -200, 5}),
    std::make_tuple("AllZeroes", std::vector<int>{0, 0, 0, 0, 0}),
    std::make_tuple("PowersOfTwo", std::vector<int>{1024, 256, 512, 128, 64, 32, 16, 8, 4, 2, 1}),
    std::make_tuple("BigNums", std::vector<int>{3243423, -1221313, 2929299, -482348, 2342453, -9876543})
};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<ChernovTRadixSortSEQ, InType>(kTestParam, PPC_SETTINGS_chernov_t_radix_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = ChernovTRadixSortFuncTests::PrintFuncTestName<ChernovTRadixSortFuncTests>;

INSTANTIATE_TEST_SUITE_P(ChernovTRadixSortTests, ChernovTRadixSortFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace chernov_t_radix_sort
