#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>

#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq/common/include/common.hpp"
#include "luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq {

class LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTestsThreads
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  // Убран static: переопределение виртуального метода
  bool CheckTestOutputData(OutType &output_data) final {
    if (!std::isfinite(output_data)) {
      return false;
    }

    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    InType test_id = std::get<0>(params);

    static const std::map<InType, double> kExpectedSums = {{2, 8.0},     {4, 64.0},    {6, 216.0},
                                                           {8, 512.0},   {10, 1000.0}, {12, 1728.0},
                                                           {16, 4096.0}, {20, 8000.0}, {24, 13824.0}};

    auto it = kExpectedSums.find(test_id);
    if (it == kExpectedSums.end()) {
      return true;
    }

    double expected = it->second;
    double rel_tol = 1e-8;
    double abs_tol = 1e-6;
    return std::abs(output_data - expected) <= ((std::abs(expected) * rel_tol) + abs_tol);
  }

  // Убран static: переопределение виртуального метода
  InType GetTestInputData() final {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    return std::get<0>(params);
  }
};

namespace {

// NOLINTNEXTLINE(readability-named-parameter)
TEST_P(LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTestsThreads, VerifyFoxMultiplication) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 9> kTestCases = {
    std::make_tuple(2, "small_even"),         std::make_tuple(4, "medium_even"),
    std::make_tuple(6, "medium_div3"),        std::make_tuple(8, "standard_block"),
    std::make_tuple(10, "large_prime_like"),  std::make_tuple(12, "large_div4"),
    std::make_tuple(16, "power_of_two"),      std::make_tuple(20, "multiple_of_block"),
    std::make_tuple(24, "max_block_coverage")};

const auto kTasks = ppc::util::AddFuncTask<LuchnikovEMultOfDenseMatrixFoxAlgoritmSeq, InType>(
    kTestCases, PPC_SETTINGS_luchnikov_e_mult_of_dense_matrices_fox_algorithm);

const auto kGtestValues = ppc::util::ExpandToValues(kTasks);
const auto kPrinter = LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTestsThreads::PrintFuncTestName<
    LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(FoxAlgorithmValidation, LuchnikovEMultOfDenseMatrixFoxAlgoritmFuncTestsThreads, kGtestValues,
                         kPrinter);

}  // namespace
}  // namespace luchnikov_e_mult_of_dense_matrices_fox_algorithm_seq
