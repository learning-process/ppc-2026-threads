#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "sabutay_sparse_complex_ccs_mult/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sabutay_sparse_complex_ccs_mult {

namespace {

// Helper function to create a random sparse matrix
CCS CreateRandomSparseMatrix(int rows, int cols, double density = 0.1) {
  CCS matrix;
  matrix.m = rows;
  matrix.n = cols;

  // Initialize col_ptr with zeros
  matrix.col_ptr.assign(cols + 1, 0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> value_dist(-10.0, 10.0);
  std::uniform_int_distribution<int> row_dist(0, rows - 1);

  int total_elements = rows * cols * density;

  for (int col = 0; col < cols; ++col) {
    int elements_in_col = (total_elements * (col + 1.0) / cols) - (total_elements * col / cols);

    for (int i = 0; i < elements_in_col; ++i) {
      int row = row_dist(gen);
      double real_part = value_dist(gen);
      double imag_part = value_dist(gen);

      matrix.row_ind.push_back(row);
      matrix.values.emplace_back(real_part, imag_part);
    }

    matrix.col_ptr[col + 1] = matrix.row_ind.size();
  }

  return matrix;
}

}  // namespace

class SabutayARunPerfTestsSeq : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    // Create test matrices with appropriate sizes for performance testing
    matrix_a_ = CreateRandomSparseMatrix(100, 100, 0.05);  // 5% density
    matrix_b_ = CreateRandomSparseMatrix(100, 100, 0.05);  // 5% density
    input_data_ = std::make_tuple(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // For performance tests, we just need to verify the output has the correct dimensions
    const CCS &a = std::get<0>(input_data_);
    const CCS &b = std::get<1>(input_data_);

    // Check that the result matrix has the correct dimensions
    return (output_data.m == a.m) && (output_data.n == b.n);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  CCS matrix_a_;
  CCS matrix_b_;
  InType input_data_;
};

TEST_P(SabutayARunPerfTestsSeq, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SabutaySparseComplexCcsMultSEQ>(PPC_SETTINGS_sabutay_sparse_complex_ccs_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SabutayARunPerfTestsSeq::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SabutayARunPerfTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult
