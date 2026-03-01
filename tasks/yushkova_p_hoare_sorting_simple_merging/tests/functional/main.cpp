#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/common/include/common.hpp"
#include "yushkova_p_hoare_sorting_simple_merging/seq/include/ops_seq.hpp"

namespace yushkova_p_hoare_sorting_simple_merging {

class YushkovaPRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param) + "_n" + std::to_string(std::get<0>(test_param).size());
  }

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }
    if (!std::ranges::is_sorted(output_data)) {
      return false;
    }

    std::vector<int> expected = input_data_;
    std::ranges::sort(expected);
    return output_data == expected;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
};

namespace {

TEST_P(YushkovaPRunFuncTestsThreads, HoareSortSimpleMergingSEQ) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    std::make_tuple(std::vector<int>{42}, "single"),
    std::make_tuple(std::vector<int>{2, 1}, "two_elements"),
    std::make_tuple(std::vector<int>{1, 2, 3, 4, 5}, "already_sorted"),
    std::make_tuple(std::vector<int>{5, 4, 3, 2, 1}, "reverse_sorted"),
    std::make_tuple(std::vector<int>{5, 1, 5, 1, 5, 1}, "duplicates"),
    std::make_tuple(std::vector<int>{0, -1, 7, -5, 2, -3}, "mixed_signs"),
    std::make_tuple(std::vector<int>{10, 3, 8, 6, 4, 9, 2, 7, 1, 5}, "random_10"),
    std::make_tuple(std::vector<int>{100, 1, 50, 2, 75, 3, 60, 4, 20, 5, 30}, "odd_count"),
    std::make_tuple(std::vector<int>{9, 9, 8, 8, 7, 7, 6, 6, 5, 5}, "pair_duplicates"),
    std::make_tuple(std::vector<int>{1000, -1000, 500, -500, 0, 250, -250}, "wide_range")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<YushkovaPHoareSortingSimpleMergingSEQ, InType>(
        kTestParam, PPC_SETTINGS_yushkova_p_hoare_sorting_simple_merging);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = YushkovaPRunFuncTestsThreads::PrintFuncTestName<YushkovaPRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(HoareSortSimpleMergingTests, YushkovaPRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace yushkova_p_hoare_sorting_simple_merging
