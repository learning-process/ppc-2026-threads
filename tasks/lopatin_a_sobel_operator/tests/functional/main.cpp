#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "lopatin_a_sobel_operator/common/include/common.hpp"
#include "lopatin_a_sobel_operator/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace lopatin_a_sobel_operator {

class LopatinARunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {};
};

namespace {

TEST_P(LopatinARunFuncTests, SobelOperator) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {"1", "2", "3"};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<LopatinASobelOperatorSEQ, InType>(kTestParam, PPC_SETTINGS_lopatin_a_sobel_operator));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = LopatinARunFuncTests::PrintFuncTestName<LopatinARunFuncTests>;

INSTANTIATE_TEST_SUITE_P(SobelOperatorTests, LopatinARunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace lopatin_a_sobel_operator
