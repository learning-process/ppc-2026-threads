#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "gasenin_l_djstra_tbb/common/include/common.hpp"
#include "gasenin_l_djstra_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace gasenin_l_djstra {

class GaseninLDjstraTbbFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

TEST_P(GaseninLDjstraTbbFuncTests, DijkstraFromParams) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<GaseninLDjstraTBB, InType>(kTestParam, PPC_SETTINGS_gasenin_l_djstra_tbb));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = GaseninLDjstraTbbFuncTests::PrintFuncTestName<GaseninLDjstraTbbFuncTests>;

INSTANTIATE_TEST_SUITE_P(DijkstraTbbTests, GaseninLDjstraTbbFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gasenin_l_djstra
