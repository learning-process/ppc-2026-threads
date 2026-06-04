#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <utility>
#include <vector>

#include "performance/include/performance.hpp"
#include "safaryan_a_sparse_matrix_mult_crs/all/include/ops_all.hpp"
#include "safaryan_a_sparse_matrix_mult_crs/common/include/common.hpp"
#include "safaryan_a_sparse_matrix_mult_crs/omp/include/ops_omp.hpp"
#include "safaryan_a_sparse_matrix_mult_crs/seq/include/ops_seq.hpp"
#include "safaryan_a_sparse_matrix_mult_crs/stl/include/ops_stl.hpp"
#include "safaryan_a_sparse_matrix_mult_crs/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace safaryan_a_sparse_matrix_mult_crs {

namespace {
SparseMatrixCCS CreateTestMatrix(int size, double density) {
  SparseMatrixCCS matrix(size, size);
  std::random_device rd;
  std::mt19937 gen(rd());
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

std::vector<std::vector<double>> DenseMultiply(const std::vector<std::vector<double>> &a,
                                               const std::vector<std::vector<double>> &b) {
  const int m = static_cast<int>(a.size());
  const int n = static_cast<int>(b[0].size());
  const int k = static_cast<int>(b.size());
  std::vector<std::vector<double>> result(m, std::vector<double>(n, 0.0));

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int idx = 0; idx < k; ++idx) {
        result[i][j] += a[i][idx] * b[idx][j];
      }
    }
  }
  return result;
}

std::vector<std::vector<double>> SparseToDense(const SparseMatrixCCS &matrix) {
  std::vector<std::vector<double>> dense(matrix.rows, std::vector<double>(matrix.cols, 0.0));
  for (int j = 0; j < matrix.cols; ++j) {
    for (int idx = matrix.col_ptrs[j]; idx < matrix.col_ptrs[j + 1]; ++idx) {
      dense[matrix.row_indices[idx]][j] = matrix.values[idx];
    }
  }
  return dense;
}
}  // namespace

class SafaryanARunPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kMatrixSize = 400;
  static constexpr double kDensity = 0.1;

  InType input_data_;
  std::vector<std::vector<double>> expected_dense_result_;

  void SetUp() override {
    const SparseMatrixCCS a = CreateTestMatrix(kMatrixSize, kDensity);
    const SparseMatrixCCS b = CreateTestMatrix(kMatrixSize, kDensity);
    input_data_ = std::make_pair(a, b);
    expected_dense_result_ = DenseMultiply(SparseToDense(a), SparseToDense(b));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != kMatrixSize || output_data.cols != kMatrixSize) {
      return false;
    }
    const std::vector<std::vector<double>> dense_result = SparseToDense(output_data);
    const double epsilon = 1e-8;
    for (int i = 0; i < kMatrixSize; ++i) {
      for (int j = 0; j < kMatrixSize; ++j) {
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

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) final {
    perf_attrs.num_running = 5;
    const auto start = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [start] {
      const auto now = std::chrono::high_resolution_clock::now();
      const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start).count();
      return static_cast<double>(ns) * 1e-9;
    };
  }
};

TEST_P(SafaryanARunPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SafaryanATaskALL, SafaryanATaskOMP, SafaryanATaskSEQ, SafaryanATaskSTL,
                                SafaryanATaskTBB>(PPC_SETTINGS_safaryan_a_sparse_matrix_mult_crs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = SafaryanARunPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultPerfTests, SafaryanARunPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace safaryan_a_sparse_matrix_mult_crs
