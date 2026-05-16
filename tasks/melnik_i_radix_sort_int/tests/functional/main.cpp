#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "melnik_i_radix_sort_int/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace melnik_i_radix_sort_int {

class MelnikIRadixSortIntFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected = input_data_;
    std::ranges::sort(expected);
    return expected == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(MelnikIRadixSortIntFuncTests, SortMatchesStdSort) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(std::vector<int>{1, 2, 3, 4, 5}, std::string("sorted_positive")),
    std::make_tuple(std::vector<int>{5, 4, 3, 2, 1}, std::string("reverse_sorted")),
    std::make_tuple(std::vector<int>{3, 1, 4, 1, 5, 9, 2, 6}, std::string("with_duplicates")),
    std::make_tuple(std::vector<int>{-5, -3, -1, 0, 1, 3, 5}, std::string("mixed_negative")),
    std::make_tuple(std::vector<int>{-10, -20, -30}, std::string("only_negative")),
    std::make_tuple(std::vector<int>{1, 2, 3}, std::string("small_array")),
    std::make_tuple(std::vector<int>{100, 50, 25, 75, 10}, std::string("random_small")),
};

const auto kTestTasksList =
    ppc::util::AddFuncTask<MelnikIRadixSortIntSEQ, InType>(kTestParam, PPC_SETTINGS_melnik_i_radix_sort_int);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = MelnikIRadixSortIntFuncTests::PrintFuncTestName<MelnikIRadixSortIntFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, MelnikIRadixSortIntFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace melnik_i_radix_sort_int
