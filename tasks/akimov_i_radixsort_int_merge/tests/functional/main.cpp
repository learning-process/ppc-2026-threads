#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "akimov_i_radixsort_int_merge/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace akimov_i_radixsort_int_merge {

class AkimovIRadixSortIntMergeFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(-1000, 1000);
    input_data_.resize(size);
    for (int &val : input_data_) {
      val = dist(gen);
    }

    expected_sorted_ = input_data_;
    std::ranges::sort(expected_sorted_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_sorted_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  InType expected_sorted_;
};

namespace {

TEST_P(AkimovIRadixSortIntMergeFuncTests, RadixSort) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(10, "10"), std::make_tuple(100, "100"),
                                            std::make_tuple(1000, "1000")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<AkimovIRadixSortIntMergeSEQ, InType>(kTestParam, PPC_SETTINGS_akimov_i_radixsort_int_merge));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = AkimovIRadixSortIntMergeFuncTests::PrintFuncTestName<AkimovIRadixSortIntMergeFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixSortTests, AkimovIRadixSortIntMergeFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace akimov_i_radixsort_int_merge
