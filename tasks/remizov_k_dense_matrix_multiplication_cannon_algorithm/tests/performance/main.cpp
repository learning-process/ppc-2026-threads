#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"
#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

class RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    int matrix_dim = 1024;
    int block_size = 64;

    std::vector<std::vector<double>> mat_a(matrix_dim, std::vector<double>(matrix_dim, 1.5));
    std::vector<std::vector<double>> mat_b(matrix_dim, std::vector<double>(matrix_dim, 2.0));

    input_data_ = std::make_tuple(block_size, mat_a, mat_b);

    expected_result_ = std::vector<std::vector<double>>(matrix_dim, std::vector<double>(matrix_dim, 3072.0));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_result_.size() != output_data.size()) {
      return false;
    }

    if (expected_result_.empty()) {
      return true;
    }

    if (expected_result_[0].size() != output_data[0].size()) {
      return false;
    }

    double epsilon = 1e-8;
    for (size_t i = 0; i < expected_result_.size(); ++i) {
      for (size_t j = 0; j < expected_result_[0].size(); ++j) {
        if (std::abs(expected_result_[i][j] - output_data[i][j]) > epsilon) {
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
  OutType expected_result_;
};

TEST_P(RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests, MultiplicationMatrixBlockSchemeCannonPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RemizovKDenseMatrixMultiplicationCannonAlgorithm>(
    PPC_SETTINGS_remizov_k_dense_matrix_multiplication_cannon_algorithm);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerfTests, RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
