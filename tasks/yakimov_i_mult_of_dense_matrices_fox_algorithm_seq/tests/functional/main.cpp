#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/common/include/common.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/seq/include/ops_seq.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq {

class YakimovIMultDenseFoxFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    if (!std::isfinite(output_data)) {
      return false;
    }

    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    InType test_id = std::get<0>(params);

    static const std::unordered_map<InType, std::pair<double, double>> kExpectedResults = {
        {1, {52.9498, 1e-4}}, {2, {-24.5668, 1e-4}}, {3, {-1.68714, 1e-4}}, {4, {0.0, 1e-9}}, {5, {60.678, 1e-4}}};

    auto it = kExpectedResults.find(test_id);
    if (it == kExpectedResults.end()) {
      return true;
    }

    double expected = it->second.first;
    double tolerance = it->second.second;

    return std::abs(output_data - expected) <= tolerance;
  }

  InType GetTestInputData() final {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    return std::get<0>(params);
  }
};

namespace {

TEST_P(YakimovIMultDenseFoxFuncTests, DenseMatrixMultiplicationFox) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kAllTestParam = {std::make_tuple(1, "small_2x3"), std::make_tuple(2, "medium_3x3"),
                                               std::make_tuple(3, "large_4x4"), std::make_tuple(4, "zero_matrices"),
                                               std::make_tuple(5, "negative_values")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<YakimovIMultOfDenseMatricesFoxAlgorithmSEQ, InType>(
    kAllTestParam, PPC_SETTINGS_yakimov_i_mult_of_dense_matrices_fox_algorithm_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = YakimovIMultDenseFoxFuncTests::PrintFuncTestName<YakimovIMultDenseFoxFuncTests>;

INSTANTIATE_TEST_SUITE_P(DenseMatrixMultiplicationTests, YakimovIMultDenseFoxFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq
