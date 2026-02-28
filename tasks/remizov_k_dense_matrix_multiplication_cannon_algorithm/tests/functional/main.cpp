#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"
#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

class RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::make_tuple(std::get<1>(params), std::get<2>(params), std::get<3>(params));
    res_ = std::get<4>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if ((res_.size() * res_[0].size()) != (output_data.size() * output_data[0].size())) {
      return false;
    }
    for (size_t i = 0; i < res_.size(); i++) {
      for (size_t j = 0; j < res_[0].size(); j++) {
        if (std::abs(res_[i][j] - output_data[i][j]) > 1e-10) {
          return false;
        }
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType res_;
};

namespace {

TEST_P(RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests, MultiplicationMatrixBlockSchemeCannon) {
  ExecuteTest(GetParam());
}

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RemizovKDenseMatrixMultiplicationCannonAlgorithm, InType>(
        kTestParam, PPC_SETTINGS_remizov_k_dense_matrix_multiplication_cannon_algorithm));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests::PrintFuncTestName<
    RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests>;

INSTANTIATE_TEST_SUITE_P(MultiplicationMatrixBlockSchemeCannonTests,
                         RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
