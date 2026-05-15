#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/common/include/test_generator.hpp"
#include "mityaeva_radix/omp/include/ops_omp.hpp"
#include "mityaeva_radix/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace mityaeva_radix {

class MityaevaRadixFunc : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

 protected:
  void SetUp() override {
    auto length = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = GenerateTest(length, length);
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

namespace {

TEST_P(MityaevaRadixFunc, RadixSortFunc) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {1, 2, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<MityaevaRadixSeq, InType>(kTestParam, PPC_SETTINGS_mityaeva_radix),
                   ppc::util::AddFuncTask<MityaevaRadixOmp, InType>(kTestParam, PPC_SETTINGS_mityaeva_radix));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = MityaevaRadixFunc::PrintFuncTestName<MityaevaRadixFunc>;

INSTANTIATE_TEST_SUITE_P(Mityaeva, MityaevaRadixFunc, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace mityaeva_radix
