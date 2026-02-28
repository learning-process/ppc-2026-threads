#include <gtest/gtest.h>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"
#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

class ChyokotovARunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};
  OutType expected_output_{};

  void SetUp() override {
    std::vector<double> &a = input_data_.first;
    std::vector<double> &b = input_data_.second;
    const int n = 700;
    a.resize(n * n);
    b.resize(n * n);
    expected_output_.resize(n * n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        a[(i * n) + j] = i + j;
        b[(i * n) + j] = i - j;
      }
    }
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        double sum = 0.0;
        for (size_t k = 0; k < n; k++) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        expected_output_[(i * n) + j] = sum;
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (expected_output_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ChyokotovARunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ChyokotovADenseMatMulFoxAlgorithmSEQ>(PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChyokotovARunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ChyokotovARunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
