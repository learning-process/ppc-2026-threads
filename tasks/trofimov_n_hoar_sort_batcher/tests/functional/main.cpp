#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "trofimov_n_hoar_sort_batcher/common/include/common.hpp"
#include "trofimov_n_hoar_sort_batcher/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace trofimov_n_hoar_sort_batcher {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<InType, std::string>;

class TrofimovNHoarSortBatcherFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    input_data_ = std::get<0>(std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(TrofimovNHoarSortBatcherFuncTests, RunFuncTests) {
  ExecuteTest(GetParam());
}

namespace {
const std::array<TestType, 12> kTestParam = {
    std::make_tuple(std::vector<int>{}, "empty"),
    std::make_tuple(std::vector<int>{1}, "single_element"),
    std::make_tuple(std::vector<int>{1, 2, 3, 4, 5}, "already_sorted"),
    std::make_tuple(std::vector<int>{5, 4, 3, 2, 1}, "reverse_sorted"),
    std::make_tuple(std::vector<int>{3, 3, 3, 3, 3}, "all_equal"),
    std::make_tuple(std::vector<int>{-5, -1, -3, 2, 0}, "negative_numbers"),
    std::make_tuple(std::vector<int>{10, 1, 9, 2, 8, 3, 7, 4, 6, 5}, "zigzag"),
    std::make_tuple(std::vector<int>{1000, -1000, 500, -500, 0}, "mixed_large_values"),
    std::make_tuple(std::vector<int>{1, 2, 3, 5, 4, 6, 7}, "almost_sorted"),
    std::make_tuple(std::vector<int>{42, 13, 7, 99, 0, -7, 88, 15}, "random_small"),
    std::make_tuple(std::vector<int>{2, 1}, "two_elements"),
    std::make_tuple(std::vector<int>{9, 1, 8, 2, 7, 3, 6, 4, 5}, "odd_count")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<TrofimovNHoarSortBatcherSEQ, InType>(kTestParam, PPC_SETTINGS_trofimov_n_hoar_sort_batcher);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = TrofimovNHoarSortBatcherFuncTests::PrintFuncTestName<TrofimovNHoarSortBatcherFuncTests>;

INSTANTIATE_TEST_SUITE_P(TrofimovNHoarSortBatcherTests, TrofimovNHoarSortBatcherFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace trofimov_n_hoar_sort_batcher
