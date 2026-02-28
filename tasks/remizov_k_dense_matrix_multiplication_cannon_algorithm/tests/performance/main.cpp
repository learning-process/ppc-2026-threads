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
  const int kCount_ = 1024;
  InType input_data_;
  OutType res_;

  void SetUp() override {
  }
  bool CheckTestOutputData(OutType &output_data) final {
  }
  InType GetTestInputData() final {
  }
};

TEST_P(RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests, MultiplicationMatrixBlockSchemeCannonPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RemizovKDenseMatrixMultiplicationCannonAlgorithm>(
    PPC_SETTINGS_remizov_k_dense_matrix_multiplication_cannon_algorithm);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RemizovKDenseMatrixMultiplicationCannonAlgorithmPerfTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
