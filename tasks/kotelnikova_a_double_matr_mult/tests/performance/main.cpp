#include <gtest/gtest.h>

#include <random>
#include <vector>
#include <algorithm>

#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"
#include "kotelnikova_a_double_matr_mult/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kotelnikova_a_double_matr_mult {

SparseMatrixCCS CreateTestMatrix(int size, double density) {
  SparseMatrixCCS matrix(size, size);
  std::mt19937 gen(42);
  std::uniform_real_distribution<double> value_dist(-10.0, 10.0);
  std::uniform_real_distribution<double> density_dist(0.0, 1.0);
  
  std::vector<std::vector<std::pair<int, double>>> columns(size);
  
  for (int j = 0; j < size; ++j) {
    for (int i = 0; i < size; ++i) {
      if (density_dist(gen) < density) {
        columns[j].emplace_back(i, value_dist(gen));
      }
    }
    std::sort(columns[j].begin(), columns[j].end());
  }
  
  matrix.col_ptrs[0] = 0;
  for (int j = 0; j < size; ++j) {
    matrix.col_ptrs[j + 1] = matrix.col_ptrs[j] + static_cast<int>(columns[j].size());
    for (const auto &[row, val] : columns[j]) {
      matrix.row_indices.push_back(row);
      matrix.values.push_back(val);
    }
  }
  
  return matrix;
}

std::vector<std::vector<double>> DenseMultiply(const std::vector<std::vector<double>> &A, 
                                               const std::vector<std::vector<double>> &B) {
  int m = static_cast<int>(A.size());
  int n = static_cast<int>(B[0].size());
  int k = static_cast<int>(B.size());
  
  std::vector<std::vector<double>> C(m, std::vector<double>(n, 0.0));
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int t = 0; t < k; ++t) {
        C[i][j] += A[i][t] * B[t][j];
      }
    }
  }
  
  return C;
}

std::vector<std::vector<double>> SparseToDense(const SparseMatrixCCS &matrix) {
  std::vector<std::vector<double>> dense(matrix.rows, std::vector<double>(matrix.cols, 0.0));
  
  for (int j = 0; j < matrix.cols; ++j) {
    for (int idx = matrix.col_ptrs[j]; idx < matrix.col_ptrs[j + 1]; ++idx) {
      int i = matrix.row_indices[idx];
      dense[i][j] = matrix.values[idx];
    }
  }
  
  return dense;
}

class KotelnikovaARunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kMatrixSize_ = 400;
  const double kDensity_ = 0.1;
  InType input_data_{};
  std::vector<std::vector<double>> expected_dense_result_;

  void SetUp() override {
    auto A = CreateTestMatrix(kMatrixSize_, kDensity_);
    auto B = CreateTestMatrix(kMatrixSize_, kDensity_);
    input_data_ = std::make_pair(A, B);
    
    auto dense_A = SparseToDense(A);
    auto dense_B = SparseToDense(B);
    expected_dense_result_ = DenseMultiply(dense_A, dense_B);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != kMatrixSize_ || output_data.cols != kMatrixSize_) {
      return false;
    }
    
    auto dense_result = SparseToDense(output_data);
    
    const double epsilon = 1e-8;
    for (int i = 0; i < kMatrixSize_; ++i) {
      for (int j = 0; j < kMatrixSize_; ++j) {
        if (std::abs(dense_result[i][j] - expected_dense_result_[i][j]) > epsilon) {
          return false;
        }
      }
    }
    
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KotelnikovaARunPerfTestSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KotelnikovaATaskSEQ>(PPC_SETTINGS_kotelnikova_a_double_matr_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KotelnikovaARunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultPerfTests, KotelnikovaARunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kotelnikova_a_double_matr_mult