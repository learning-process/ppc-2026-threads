#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <tuple>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {
int64_t run_and_time(const CRSMatrix &a, const CRSMatrix &b, CRSMatrix *out) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));

  EXPECT_TRUE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());

  const auto start = std::chrono::high_resolution_clock::now();
  EXPECT_TRUE(task.Run());
  const auto end = std::chrono::high_resolution_clock::now();

  EXPECT_TRUE(task.PostProcessing());

  *out = task.GetOutput();
  return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

CRSMatrix build_band_matrix(int n, const Complex &value, int bandwidth = 5) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    for (int offset = 0; offset < bandwidth; ++offset) {
      int col = i + offset;
      if (col < n) {
        m.col_indices.push_back(col);
        m.values.push_back(value);
      }
    }
    m.row_ptr[i + 1] = static_cast<int>(m.col_indices.size());
  }
  return m;
}

CRSMatrix build_diagonal_matrix(int n, const Complex &value) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    m.col_indices.push_back(i);
    m.values.push_back(value);
    m.row_ptr[i + 1] = i + 1;
  }
  return m;
}

CRSMatrix build_scattered_matrix(int n, const Complex &value, int step = 7) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    int col = (i * step) % n;
    m.col_indices.push_back(col);
    m.values.push_back(value);
    m.row_ptr[i + 1] = i + 1;
  }

  return m;
}

CRSMatrix build_block_diagonal_matrix(int n, int block_size, const Complex &value) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    int block_start = (i / block_size) * block_size;
    int block_end = std::min(block_start + block_size, n);
    for (int col = block_start; col < block_end; ++col) {
      m.col_indices.push_back(col);
      m.values.push_back(value);
    }
    m.row_ptr[i + 1] = static_cast<int>(m.col_indices.size());
  }
  return m;
}
}  // namespace

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, MeasureTime) {
  const int n = 1000;

  CRSMatrix a = build_band_matrix(n, Complex(1.0, 1.0));
  CRSMatrix b = build_band_matrix(n, Complex(1.0, 0.0));

  ASSERT_TRUE(a.is_valid()) << "Matrix a failed CRS validation";
  ASSERT_TRUE(b.is_valid()) << "Matrix b failed CRS validation";

  CRSMatrix c;
  const int64_t duration_ms = run_and_time(a, b, &c);
  EXPECT_LT(duration_ms, 5000) << "Multiplication took " << duration_ms << " ms, which exceeds the 5 s limit";

  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, n);
  EXPECT_EQ(c.cols, n);
  EXPECT_GT(c.non_zeros(), 0U);
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, DiagonalMatrix5000x5000) {
  const int n = 5000;

  CRSMatrix a = build_diagonal_matrix(n, Complex(2.0, 1.0));
  CRSMatrix b = build_diagonal_matrix(n, Complex(1.0, -1.0));

  ASSERT_TRUE(a.is_valid());
  ASSERT_TRUE(b.is_valid());

  CRSMatrix c;
  const int64_t duration_ms = run_and_time(a, b, &c);
  EXPECT_LT(duration_ms, 500) << "Diagonal 5000x5000 took " << duration_ms << " ms, expected under 500 ms";

  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, n);
  EXPECT_EQ(c.cols, n);

  EXPECT_EQ(c.non_zeros(), static_cast<size_t>(n));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, WideBandMatrix2000x2000) {
  const int n = 2000;

  CRSMatrix a = build_band_matrix(n, Complex(1.0, 0.5), 20);
  CRSMatrix b = build_band_matrix(n, Complex(0.5, 1.0), 20);

  ASSERT_TRUE(a.is_valid());
  ASSERT_TRUE(b.is_valid());

  CRSMatrix c;
  const int64_t duration_ms = run_and_time(a, b, &c);
  EXPECT_LT(duration_ms, 5000) << "Wide-band 2000x2000 took " << duration_ms << " ms, exceeds 5 s limit";

  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, n);
  EXPECT_EQ(c.cols, n);
  EXPECT_GT(c.non_zeros(), 0U);
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, ScatteredMatrix3000x3000) {
  const int n = 3000;

  CRSMatrix a = build_scattered_matrix(n, Complex(1.0, 1.0), 7);
  CRSMatrix b = build_scattered_matrix(n, Complex(1.0, -1.0), 11);

  ASSERT_TRUE(a.is_valid());
  ASSERT_TRUE(b.is_valid());

  CRSMatrix c;
  const int64_t duration_ms = run_and_time(a, b, &c);
  EXPECT_LT(duration_ms, 1000) << "Scattered 3000x3000 took " << duration_ms << " ms, expected under 1 s";

  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, n);
  EXPECT_EQ(c.cols, n);
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, BlockDiagonalMatrix1000x1000) {
  const int n = 1000;
  const int block_size = 10;

  CRSMatrix a = build_block_diagonal_matrix(n, block_size, Complex(1.0, 1.0));
  CRSMatrix b = build_block_diagonal_matrix(n, block_size, Complex(2.0, -1.0));

  ASSERT_TRUE(a.is_valid());
  ASSERT_TRUE(b.is_valid());

  CRSMatrix c;
  const int64_t duration_ms = run_and_time(a, b, &c);
  EXPECT_LT(duration_ms, 2000) << "Block-diagonal 1000x1000 took " << duration_ms << " ms, exceeds 2 s limit";

  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, n);
  EXPECT_EQ(c.cols, n);

  EXPECT_EQ(c.non_zeros(), static_cast<size_t>(n * block_size));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, RectangularBandMatrix500x2000x500) {
  const int m = 500;
  const int k = 2000;
  const int n_out = 500;

  CRSMatrix a(m, k);
  for (int i = 0; i < m; ++i) {
    for (int offset = 0; offset < 5; ++offset) {
      int col = (i * (k / m)) + offset;
      if (col < k) {
        a.col_indices.push_back(col);
        a.values.emplace_back(1.0, 0.5);
      }
    }
    a.row_ptr[i + 1] = static_cast<int>(a.col_indices.size());
  }

  CRSMatrix b(k, n_out);
  for (int i = 0; i < k; ++i) {
    int col = i * n_out / k;
    if (col < n_out) {
      b.col_indices.push_back(col);
      b.values.emplace_back(0.5, 1.0);
    }
    b.row_ptr[i + 1] = static_cast<int>(b.col_indices.size());
  }

  ASSERT_TRUE(a.is_valid());
  ASSERT_TRUE(b.is_valid());

  CRSMatrix c;
  const int64_t duration_ms = run_and_time(a, b, &c);
  EXPECT_LT(duration_ms, 2000) << "Rectangular 500x2000x500 took " << duration_ms << " ms, exceeds 2 s limit";

  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, m);
  EXPECT_EQ(c.cols, n_out);
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
