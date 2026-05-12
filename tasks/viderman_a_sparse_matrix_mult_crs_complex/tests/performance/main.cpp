#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <tuple>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {
struct RunResult {
  bool ok = false;
  int64_t duration_ms = 0;
  CRSMatrix out;
};

RunResult RunAndTime(const CRSMatrix &a, const CRSMatrix &b) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  if (!task.Validation() || !task.PreProcessing()) {
    return {};
  }

  const auto start = std::chrono::high_resolution_clock::now();
  const bool run_ok = task.Run();
  const auto end = std::chrono::high_resolution_clock::now();

  if (!run_ok || !task.PostProcessing()) {
    return {};
  }

  RunResult result;
  result.ok = true;
  result.duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  result.out = task.GetOutput();
  return result;
}

bool CheckResult(const RunResult &result, int64_t limit_ms, int rows, int cols, std::size_t min_nnz,
                 std::optional<std::size_t> exact_nnz = std::nullopt) {
  if (!result.ok) {
    return false;
  }
  if (result.duration_ms >= limit_ms) {
    return false;
  }
  if (!result.out.IsValid() || result.out.rows != rows || result.out.cols != cols) {
    return false;
  }
  if (result.out.NonZeros() < min_nnz) {
    return false;
  }
  if (exact_nnz.has_value() && result.out.NonZeros() != exact_nnz.value()) {
    return false;
  }
  return true;
}

CRSMatrix BuildBandMatrix(int n, const Complex &value, int bandwidth = 5) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    for (int offset = 0; offset < bandwidth; ++offset) {
      const int col = i + offset;
      if (col < n) {
        m.col_indices.push_back(col);
        m.values.push_back(value);
      }
    }
    m.row_ptr[i + 1] = static_cast<int>(m.col_indices.size());
  }
  return m;
}

CRSMatrix BuildDiagonalMatrix(int n, const Complex &value) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    m.col_indices.push_back(i);
    m.values.push_back(value);
    m.row_ptr[i + 1] = i + 1;
  }
  return m;
}

CRSMatrix BuildScatteredMatrix(int n, const Complex &value, int step = 7) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const int col = (i * step) % n;
    m.col_indices.push_back(col);
    m.values.push_back(value);
    m.row_ptr[i + 1] = i + 1;
  }
  return m;
}

CRSMatrix BuildBlockDiagonalMatrix(int n, int block_size, const Complex &value) {
  CRSMatrix m(n, n);
  for (int i = 0; i < n; ++i) {
    const int block_start = (i / block_size) * block_size;
    const int block_end = std::min(block_start + block_size, n);
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
  const RunResult result = RunAndTime(BuildBandMatrix(n, Complex(1.0, 1.0)), BuildBandMatrix(n, Complex(1.0, 0.0)));
  EXPECT_TRUE(CheckResult(result, 5000, n, n, 1U));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, DiagonalMatrix5000x5000) {
  const int n = 5000;
  const RunResult result =
      RunAndTime(BuildDiagonalMatrix(n, Complex(2.0, 1.0)), BuildDiagonalMatrix(n, Complex(1.0, -1.0)));
  EXPECT_TRUE(CheckResult(result, 500, n, n, static_cast<std::size_t>(n), static_cast<std::size_t>(n)));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, WideBandMatrix2000x2000) {
  const int n = 2000;
  const RunResult result =
      RunAndTime(BuildBandMatrix(n, Complex(1.0, 0.5), 20), BuildBandMatrix(n, Complex(0.5, 1.0), 20));
  EXPECT_TRUE(CheckResult(result, 5000, n, n, 1U));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, ScatteredMatrix3000x3000) {
  const int n = 3000;
  const RunResult result =
      RunAndTime(BuildScatteredMatrix(n, Complex(1.0, 1.0), 7), BuildScatteredMatrix(n, Complex(1.0, -1.0), 11));
  EXPECT_TRUE(CheckResult(result, 1000, n, n, 1U));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, BlockDiagonalMatrix1000x1000) {
  const int n = 1000;
  const int block_size = 10;
  const RunResult result = RunAndTime(BuildBlockDiagonalMatrix(n, block_size, Complex(1.0, 1.0)),
                                      BuildBlockDiagonalMatrix(n, block_size, Complex(2.0, -1.0)));
  EXPECT_TRUE(CheckResult(result, 2000, n, n, static_cast<std::size_t>(n * block_size),
                          static_cast<std::size_t>(n * block_size)));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, RectangularBandMatrix500x2000x500) {
  const int m = 500;
  const int k = 2000;
  const int n_out = 500;

  CRSMatrix a(m, k);
  for (int i = 0; i < m; ++i) {
    for (int offset = 0; offset < 5; ++offset) {
      const int col = (i * (k / m)) + offset;
      if (col < k) {
        a.col_indices.push_back(col);
        a.values.emplace_back(1.0, 0.5);
      }
    }
    a.row_ptr[i + 1] = static_cast<int>(a.col_indices.size());
  }

  CRSMatrix b(k, n_out);
  for (int i = 0; i < k; ++i) {
    const int col = (i * n_out) / k;
    if (col < n_out) {
      b.col_indices.push_back(col);
      b.values.emplace_back(0.5, 1.0);
    }
    b.row_ptr[i + 1] = static_cast<int>(b.col_indices.size());
  }

  const RunResult result = RunAndTime(a, b);
  EXPECT_TRUE(CheckResult(result, 2000, m, n_out, 1U));
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
