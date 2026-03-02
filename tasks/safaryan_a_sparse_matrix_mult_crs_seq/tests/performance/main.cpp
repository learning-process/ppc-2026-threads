#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "safaryan_a_sparse_matrix_mult_crs_seq/common/include/common.hpp"
#include "safaryan_a_sparse_matrix_mult_crs_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {

namespace {

CRSMatrix CreateTestMatrix(size_t size, double density) {
  CRSMatrix matrix;
  matrix.rows = size;
  matrix.cols = size;
  matrix.row_ptr.assign(size + 1, 0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> value_dist(-10.0, 10.0);
  std::uniform_real_distribution<double> density_dist(0.0, 1.0);

  // build per-row list of (col, val)
  std::vector<std::vector<std::pair<size_t, double>>> rows(size);

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      if (density_dist(gen) < density) {
        rows[i].emplace_back(j, value_dist(gen));
      }
    }
    std::sort(rows[i].begin(), rows[i].end(), [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
  }

  matrix.row_ptr[0] = 0;
  for (size_t i = 0; i < size; ++i) {
    for (const auto &[col, val] : rows[i]) {
      matrix.col_indices.push_back(col);
      matrix.values.push_back(val);
    }
    matrix.row_ptr[i + 1] = matrix.values.size();
  }

  matrix.nnz = matrix.values.size();
  return matrix;
}

std::vector<std::vector<double>> DenseMultiply(const std::vector<std::vector<double>> &a,
                                               const std::vector<std::vector<double>> &b) {
  const size_t m = a.size();
  const size_t n = b.empty() ? 0 : b[0].size();
  const size_t k = b.size();

  std::vector<std::vector<double>> result(m, std::vector<double>(n, 0.0));

  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      for (size_t idx = 0; idx < k; ++idx) {
        result[i][j] += a[i][idx] * b[idx][j];
      }
    }
  }

  return result;
}

std::vector<std::vector<double>> SparseToDense(const CRSMatrix &matrix) {
  std::vector<std::vector<double>> dense(matrix.rows, std::vector<double>(matrix.cols, 0.0));

  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t idx = matrix.row_ptr[i]; idx < matrix.row_ptr[i + 1]; ++idx) {
      dense[i][matrix.col_indices[idx]] = matrix.values[idx];
    }
  }

  return dense;
}

}  // namespace

class SafaryanARunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kMatrixSize = 400;
  static constexpr double kDensity = 0.1;

  InType input_data_;
  std::vector<std::vector<double>> expected_dense_result_;

  void SetUp() override {
    const CRSMatrix a = CreateTestMatrix(kMatrixSize, kDensity);
    const CRSMatrix b = CreateTestMatrix(kMatrixSize, kDensity);
    input_data_ = std::make_pair(a, b);

    const std::vector<std::vector<double>> dense_a = SparseToDense(a);
    const std::vector<std::vector<double>> dense_b = SparseToDense(b);
    expected_dense_result_ = DenseMultiply(dense_a, dense_b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != kMatrixSize || output_data.cols != kMatrixSize) {
      return false;
    }

    const std::vector<std::vector<double>> dense_result = SparseToDense(output_data);

    const double epsilon = 1e-8;
    for (size_t i = 0; i < kMatrixSize; ++i) {
      for (size_t j = 0; j < kMatrixSize; ++j) {
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

TEST_P(SafaryanARunPerfTestSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SafaryanATaskSEQ>(PPC_SETTINGS_safaryan_a_sparse_matrix_mult_crs_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SafaryanARunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultPerfTests, SafaryanARunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
