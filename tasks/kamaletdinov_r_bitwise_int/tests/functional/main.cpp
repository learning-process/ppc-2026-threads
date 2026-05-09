#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "kamaletdinov_r_bitwise_int/all/include/ops_all.hpp"
#include "kamaletdinov_r_bitwise_int/common/include/common.hpp"
#include "kamaletdinov_r_bitwise_int/omp/include/ops_omp.hpp"
#include "kamaletdinov_r_bitwise_int/seq/include/ops_seq.hpp"
#include "kamaletdinov_r_bitwise_int/stl/include/ops_stl.hpp"
#include "kamaletdinov_r_bitwise_int/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kamaletdinov_r_bitwise_int {

class KamaletdinovRBitwiseIntRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(KamaletdinovRBitwiseIntRunFuncTests, BitwiseSort) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {std::make_tuple(0, "empty"),      std::make_tuple(1, "single"),
                                             std::make_tuple(2, "two"),        std::make_tuple(3, "three"),
                                             std::make_tuple(5, "five"),       std::make_tuple(7, "seven"),
                                             std::make_tuple(10, "ten"),       std::make_tuple(100, "hundred"),
                                             std::make_tuple(256, "pow2_256"), std::make_tuple(1000, "thousand")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KamaletdinovRBitwiseIntALL, InType>(kTestParam, PPC_SETTINGS_kamaletdinov_r_bitwise_int),
    ppc::util::AddFuncTask<KamaletdinovRBitwiseIntOMP, InType>(kTestParam, PPC_SETTINGS_kamaletdinov_r_bitwise_int),
    ppc::util::AddFuncTask<KamaletdinovRBitwiseIntSEQ, InType>(kTestParam, PPC_SETTINGS_kamaletdinov_r_bitwise_int),
    ppc::util::AddFuncTask<KamaletdinovRBitwiseIntSTL, InType>(kTestParam, PPC_SETTINGS_kamaletdinov_r_bitwise_int),
    ppc::util::AddFuncTask<KamaletdinovRBitwiseIntTBB, InType>(kTestParam, PPC_SETTINGS_kamaletdinov_r_bitwise_int));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KamaletdinovRBitwiseIntRunFuncTests::PrintFuncTestName<KamaletdinovRBitwiseIntRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(BitwiseIntTests, KamaletdinovRBitwiseIntRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamaletdinov_r_bitwise_int
