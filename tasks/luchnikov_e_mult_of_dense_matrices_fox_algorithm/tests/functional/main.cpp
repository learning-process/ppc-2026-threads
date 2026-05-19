#include <gtest/gtest.h>

#include <array>
#include <string>
#include <tuple>

#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm/common/include/common.hpp"
#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm {

class LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
    return (output_data > 0);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTests, UnifiedMatrixMultiplication) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParams = {std::make_tuple(3, "small"), std::make_tuple(5, "medium"),
                                             std::make_tuple(7, "large")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq, InType>(
                       kTestParams, PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm),
                   ppc::util::AddFuncTask<LuchnikovEMultOfDenseMatrixFoxAlgoritmOMP, InType>(
                       kTestParams, PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm),
                   ppc::util::AddFuncTask<LuchnikovEMultOfDenseMatrixFoxAlgoritmSTL, InType>(
                       kTestParams, PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm),
                   ppc::util::AddFuncTask<LuchnikovEMultOfDenseMatrixFoxAlgoritmTBB, InType>(
                       kTestParams, PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName =
    LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTests::PrintFuncTestName<LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTests>;

INSTANTIATE_TEST_SUITE_P(AllTechTests, LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm
