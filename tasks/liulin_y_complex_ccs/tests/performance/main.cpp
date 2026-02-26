#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <random>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"
#include "liulin_y_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace liulin_y_complex_ccs {

namespace {

CCSMatrix CreateRandomSparseComplexMatrix(int rows, int cols, double density) {
  CCSMatrix matrix;
  matrix.count_rows = rows;
  matrix.count_cols = cols;
  matrix.col_index.resize(cols + 1, 0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  std::uniform_real_distribution<double> val_dist(-10.0, 10.0);

  int total_nnz = 0;
  for (int j = 0; j < cols; ++j) {
    std::vector<int> current_col_rows;
    for (int i = 0; i < rows; ++i) {
      if (prob_dist(gen) < density) {
        current_col_rows.push_back(i);
      }
    }

    for (int row : current_col_rows) {
      matrix.values.emplace_back(val_dist(gen), val_dist(gen));
      matrix.row_index.push_back(row);
      total_nnz++;
    }
    matrix.col_index[j + 1] = total_nnz;
  }

  return matrix;
}

}  // namespace

class LiulinYComplexCcsPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    int size = 1000;
    double density = 0.01;

    matrix_a_ = CreateRandomSparseComplexMatrix(size, size, density);
    matrix_b_ = CreateRandomSparseComplexMatrix(size, size, density);

    input_data_ = std::make_pair(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.count_rows != matrix_a_.count_rows || output_data.count_cols != matrix_b_.count_cols) {
      return false;
    }

    if (output_data.col_index.size() != static_cast<size_t>(output_data.count_cols + 1)) {
      return false;
    }

    for (size_t i = 0; i < output_data.col_index.size() - 1; ++i) {
      if (output_data.col_index[i] > output_data.col_index[i + 1]) {
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

TEST_P(LiulinYComplexCcsPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, LiulinYComplexCcs>(PPC_SETTINGS_liulin_y_complex_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LiulinYComplexCcsPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(MatrixMultiplyPerfTests, LiulinYComplexCcsPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace liulin_y_complex_ccs
