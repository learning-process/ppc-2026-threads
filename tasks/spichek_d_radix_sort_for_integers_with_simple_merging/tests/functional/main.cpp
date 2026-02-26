#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <random>
#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }
 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int n = std::get<0>(params);
    input_data_.resize(n);
    std::mt19937 gen(42);
    std::uniform_int_distribution<> dis(-1000, 1000);
    for(int i = 0; i < n; ++i) input_data_[i] = dis(gen);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
  }
  InType GetTestInputData() final { return input_data_; }
 private:
  InType input_data_;
};

TEST_P(RadixSortFuncTests, SequentialSort) { ExecuteTest(GetParam()); }

const std::array<TestType, 3> kTestParam = {
    std::make_tuple(10, "Small_10"),
    std::make_tuple(100, "Medium_100"),
    std::make_tuple(1000, "Large_1000")
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<RadixSortSimpleMergingSEQ, InType>(kTestParam, "spichek_d_radix_sort_for_integers_with_simple_merging")
);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
INSTANTIATE_TEST_SUITE_P(RadixSortTests, RadixSortFuncTests, kGtestValues, RadixSortFuncTests::PrintFuncTestName<RadixSortFuncTests>);

} // namespace spichek_d_radix_sort_for_integers_with_simple_merging