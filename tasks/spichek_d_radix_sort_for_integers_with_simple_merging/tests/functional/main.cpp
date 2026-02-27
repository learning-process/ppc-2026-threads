#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <ranges>
#include <string>
#include <tuple>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class SpichekDRadixSortRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int vector_size = std::get<0>(params);

    input_data.resize(vector_size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-10000, 10000);

    for (int i = 0; i < vector_size; ++i) {
      input_data[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data;
  }

  InType input_data;
};

namespace {

TEST_P(SpichekDRadixSortRunFuncTests, SpichekDRadixSortSequential) {
  ExecuteTest(GetParam());
}

TEST_P(SpichekDRadixSortRunFuncTests, SpichekDRadixSortWithNegativeValues) {
  for (auto &val : input_data) {
    val = -std::abs(val);
  }
  ExecuteTest(GetParam());
}

TEST_P(SpichekDRadixSortRunFuncTests, SpichekDRadixSortWithAllSameValues) {
  std::ranges::fill(input_data, 42);
  ExecuteTest(GetParam());
}

TEST_P(SpichekDRadixSortRunFuncTests, SpichekDRadixSortWithSortedInput) {
  std::ranges::sort(input_data);
  ExecuteTest(GetParam());
}

TEST_P(SpichekDRadixSortRunFuncTests, SpichekDRadixSortWithReversedInput) {
  std::ranges::sort(input_data, std::greater<>{});
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(100, "small_vector"),
                                            std::make_tuple(1000, "medium_vector"),
                                            std::make_tuple(5000, "large_vector")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<SpichekDRadixSortSEQ, InType>(
    kTestParam, PPC_SETTINGS_spichek_d_radix_sort_for_integers_with_simple_merging));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SpichekDRadixSortRunFuncTests::PrintFuncTestName<SpichekDRadixSortRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(SpichekDRadixSortTests, SpichekDRadixSortRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
