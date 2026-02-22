#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace barkalova_m_mult_matrix_ccs {

namespace {

CCSMatrix CreateRandomSparseMatrix(int rows, int cols, double density) {
  CCSMatrix matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.col_ptrs.resize(cols + 1, 0);

  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  static std::uniform_real_distribution<double> value_dist(0.1, 10.0);

  std::vector<int> col_counts(cols, 0);
  std::vector<std::vector<Complex>> col_values(cols);
  std::vector<std::vector<int>> col_rows(cols);

  int total_nnz = 0;

  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      double random_value = dist(gen);

      if (random_value < density) {
        double real = value_dist(gen);
        double imag = value_dist(gen);

        col_values[col].emplace_back(real, imag);
        col_rows[col].push_back(row);
        ++col_counts[col];
        ++total_nnz;
      }
    }
  }

  matrix.nnz = total_nnz;
  matrix.values.resize(total_nnz);
  matrix.row_indices.resize(total_nnz);

  int current_index = 0;
  matrix.col_ptrs[0] = 0;

  for (int col = 0; col < cols; ++col) {
    for (int i = 0; i < col_counts[col]; ++i) {
      matrix.values[current_index] = col_values[col][i];
      matrix.row_indices[current_index] = col_rows[col][i];
      ++current_index;
    }
    matrix.col_ptrs[col + 1] = current_index;
  }

  return matrix;
}

}  // namespace

class BarkalovaMMultMatrixCcsPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    int size = 1500;
    double density = 0.01;

    matrix_a_ = CreateRandomSparseMatrix(size, size, density);
    matrix_b_ = CreateRandomSparseMatrix(size, size, density);

    input_data_ = std::make_pair(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != matrix_a_.rows || output_data.cols != matrix_b_.cols) {
      return false;
    }

    if (output_data.col_ptrs.size() != static_cast<size_t>(output_data.cols) + 1) {
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

TEST_P(BarkalovaMMultMatrixCcsPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BarkalovaMMultMatrixCcsSEQ>(PPC_SETTINGS_barkalova_m_mult_matrix_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BarkalovaMMultMatrixCcsPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(MatrixMultiplyPerfTests, BarkalovaMMultMatrixCcsPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace barkalova_m_mult_matrix_ccs
