#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "safaryan_a_sparse_matrix_mult_crs_seq/common/include/common.hpp"
#include "safaryan_a_sparse_matrix_mult_crs_seq/seq/include/ops_seq.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {
namespace {

// ------------------ helpers: Dense <-> CRS ------------------

CRSMatrix DenseToCrs(const std::vector<std::vector<double>> &dense, double eps = 0.0) {
  CRSMatrix m;
  m.rows = dense.size();
  m.cols = dense.empty() ? 0 : dense[0].size();
  m.row_ptr.assign(m.rows + 1, 0);

  for (size_t i = 0; i < m.rows; ++i) {
    m.row_ptr[i] = m.values.size();
    for (size_t j = 0; j < m.cols; ++j) {
      const double v = dense[i][j];
      if (std::abs(v) > eps) {
        m.values.push_back(v);
        m.col_indices.push_back(j);
      }
    }
  }

  m.nnz = m.values.size();
  m.row_ptr[m.rows] = m.nnz;
  return m;
}

std::vector<std::vector<double>> CrsToDense(const CRSMatrix &m) {
  std::vector<std::vector<double>> dense(m.rows, std::vector<double>(m.cols, 0.0));

  for (size_t i = 0; i < m.rows; ++i) {
    for (size_t idx = m.row_ptr[i]; idx < m.row_ptr[i + 1]; ++idx) {
      dense[i][m.col_indices[idx]] += m.values[idx];
    }
  }

  return dense;
}

std::vector<std::vector<double>> DenseMul(const std::vector<std::vector<double>> &a,
                                          const std::vector<std::vector<double>> &b) {
  const size_t n = a.size();
  const size_t k = a.empty() ? 0 : a[0].size();
  const size_t m = b.empty() ? 0 : b[0].size();

  std::vector<std::vector<double>> c(n, std::vector<double>(m, 0.0));
  for (size_t i = 0; i < n; ++i) {
    for (size_t kk = 0; kk < k; ++kk) {
      const double av = a[i][kk];
      if (av == 0.0) {
        continue;
      }
      for (size_t j = 0; j < m; ++j) {
        c[i][j] += av * b[kk][j];
      }
    }
  }

  return c;
}

void ExpectDenseEqual(const std::vector<std::vector<double>> &x, const std::vector<std::vector<double>> &y,
                      double eps = 1e-9) {
  ASSERT_EQ(x.size(), y.size());
  if (x.empty()) {
    return;
  }
  ASSERT_EQ(x[0].size(), y[0].size());

  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < x[0].size(); ++j) {
      ASSERT_NEAR(x[i][j], y[i][j], eps) << "Mismatch at (" << i << "," << j << ")";
    }
  }
}

// ------------------ random generator ------------------

std::vector<std::vector<double>> GenDense(size_t rows, size_t cols, double density, std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> prob(0.0, 1.0);
  std::uniform_real_distribution<double> val(-5.0, 5.0);

  std::vector<std::vector<double>> d(rows, std::vector<double>(cols, 0.0));
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      if (prob(rng) < density) {
        d[i][j] = val(rng);
      }
    }
  }

  return d;
}

// -------------- functional tests --------------

TEST(SafaryanASparseMatrixMultCRSSeqFunctional, SmallFixedCase) {
  const std::vector<std::vector<double>> a_dense = {
      {1.0, 0.0, 2.0},
      {0.0, -1.0, 3.0},
  };
  const std::vector<std::vector<double>> b_dense = {
      {2.0, 1.0},
      {0.0, 4.0},
      {1.0, -2.0},
  };

  const CRSMatrix a_crs = DenseToCrs(a_dense);
  const CRSMatrix b_crs = DenseToCrs(b_dense);

  SafaryanASparseMatrixMultCRSSeq task({a_crs, b_crs});
  ASSERT_TRUE(task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing());

  const auto got = CrsToDense(task.GetOutput());
  const auto ref = DenseMul(a_dense, b_dense);

  ExpectDenseEqual(got, ref);
}

TEST(SafaryanASparseMatrixMultCRSSeqFunctional, RandomCaseMatchesDense) {
  const size_t n = 25;
  const size_t k = 30;
  const size_t m = 20;

  const auto a_dense = GenDense(n, k, 0.12, 100U);
  const auto b_dense = GenDense(k, m, 0.12, 200U);

  const CRSMatrix a_crs = DenseToCrs(a_dense);
  const CRSMatrix b_crs = DenseToCrs(b_dense);

  SafaryanASparseMatrixMultCRSSeq task({a_crs, b_crs});
  ASSERT_TRUE(task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing());

  const auto got = CrsToDense(task.GetOutput());
  const auto ref = DenseMul(a_dense, b_dense);

  ExpectDenseEqual(got, ref, 1e-8);
}

TEST(SafaryanASparseMatrixMultCRSSeqFunctional, ValidationFailsOnBadSizes) {
  CRSMatrix a_crs;
  a_crs.rows = 2;
  a_crs.cols = 3;
  a_crs.row_ptr = {0, 0, 0};
  a_crs.nnz = 0;

  CRSMatrix b_crs;
  b_crs.rows = 4;
  b_crs.cols = 2;
  b_crs.row_ptr = {0, 0, 0, 0, 0};
  b_crs.nnz = 0;

  SafaryanASparseMatrixMultCRSSeq task({a_crs, b_crs});
  ASSERT_FALSE(task.Validation());
}

}  // namespace

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
