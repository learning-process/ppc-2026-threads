#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs_seq/common/include/common.hpp"
#include "kapanova_s_sparse_matrix_mult_ccs_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs_seq {

namespace {

CCSMatrix CreateRandomSparseMatrix(size_t rows, size_t cols, double density) {
  CCSMatrix matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.col_ptrs.resize(cols + 1, 0);

  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  static std::uniform_real_distribution<double> val_dist(-10.0, 10.0);

  std::vector<size_t> col_counts(cols, 0);
  std::vector<std::vector<double>> col_values(cols);
  std::vector<std::vector<size_t>> col_rows(cols);

  size_t total_nnz = 0;

  for (size_t col = 0; col < cols; ++col) {
    for (size_t row = 0; row < rows; ++row) {
      if (prob_dist(gen) < density) {
        double value = val_dist(gen);
        col_values[col].push_back(value);
        col_rows[col].push_back(row);
        ++col_counts[col];
        ++total_nnz;
      }
    }
  }

  matrix.nnz = total_nnz;
  matrix.values.resize(total_nnz);
  matrix.row_indices.resize(total_nnz);

  size_t current_index = 0;
  matrix.col_ptrs[0] = 0;

  for (size_t col = 0; col < cols; ++col) {
    for (size_t i = 0; i < col_counts[col]; ++i) {
      matrix.values[current_index] = col_values[col][i];
      matrix.row_indices[current_index] = col_rows[col][i];
      ++current_index;
    }
    matrix.col_ptrs[col + 1] = current_index;
  }

  return matrix;
}

}  // namespace

class KapanovaSMatrixMultiplyPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    size_t size = 1500;
    double density = 0.01;

    matrix_a_ = CreateRandomSparseMatrix(size, size, density);
    matrix_b_ = CreateRandomSparseMatrix(size, size, density);

    input_data_ = std::make_pair(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != matrix_a_.rows || output_data.cols != matrix_b_.cols) {
      return false;
    }

    if (output_data.col_ptrs.size() != output_data.cols + 1) {
      return false;
    }

    for (size_t i = 0; i < output_data.col_ptrs.size() - 1; ++i) {
      if (output_data.col_ptrs[i] > output_data.col_ptrs[i + 1]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  CCSMatrix matrix_a_;
  CCSMatrix matrix_b_;
  InType input_data_;
};

TEST_P(KapanovaSMatrixMultiplyPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KapanovaSSparseMatrixMultCCSSeq>(
    PPC_SETTINGS_kapanova_s_sparse_matrix_mult_ccs_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KapanovaSMatrixMultiplyPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(MatrixMultiplyPerfTests, KapanovaSMatrixMultiplyPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kapanova_s_sparse_matrix_mult_ccs_seq
