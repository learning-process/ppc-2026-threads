#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "gasenin_l_djstra/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace gasenin_l_djstra {

class GaseninLRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = input_data_ * (input_data_ - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
  OutType expected_output_ = 0;
};

namespace {

TEST_P(GaseninLRunFuncTestsThreads, DijkstraFromParams) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<GaseninLDjstraSEQ, InType>(kTestParam, PPC_SETTINGS_gasenin_l_djstra));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = GaseninLRunFuncTestsThreads::PrintFuncTestName<GaseninLRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(DijkstraSeqTests, GaseninLRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gasenin_l_djstra
