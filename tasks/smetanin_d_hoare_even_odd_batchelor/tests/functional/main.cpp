#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>

#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"
#include "smetanin_d_hoare_even_odd_batchelor/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace smetanin_d_hoare_even_odd_batchelor {

class SmetaninDRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int size = std::get<0>(params);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(-10000, 10000);
    input_data_.resize(static_cast<std::size_t>(size));
    for (int &val : input_data_) {
      val = dist(rng);
    }
    expected_ = input_data_;
    std::ranges::sort(expected_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_;
};

namespace {

TEST_P(SmetaninDRunFuncTests, SortTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParam = {
    std::make_tuple(0, "empty"),   std::make_tuple(1, "single"), std::make_tuple(2, "two"),
    std::make_tuple(10, "small"),  std::make_tuple(100, "mid"),  std::make_tuple(1000, "large"),
    std::make_tuple(10000, "big"),
};

const auto kTestTasksList =
    ppc::util::AddFuncTask<SmetaninDHoarSortSEQ, InType>(kTestParam, PPC_SETTINGS_smetanin_d_hoare_even_odd_batchelor);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = SmetaninDRunFuncTests::PrintFuncTestName<SmetaninDRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(SortTests, SmetaninDRunFuncTests, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace smetanin_d_hoare_even_odd_batchelor
