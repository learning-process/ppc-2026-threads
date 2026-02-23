#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <vector>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

static CRSMatrix BuildBandMatrix(int n, const Complex &value, int bandwidth = 5) {
  CRSMatrix M(n, n);
  for (int i = 0; i < n; ++i) {
    for (int offset = 0; offset < bandwidth; ++offset) {
      int col = i + offset;
      if (col < n) {
        M.col_indices.push_back(col);
        M.values.push_back(value);
      }
    }
    M.row_ptr[i + 1] = static_cast<int>(M.col_indices.size());
  }
  return M;
}

static CRSMatrix BuildDiagonalMatrix(int n, const Complex &value) {
  CRSMatrix M(n, n);
  for (int i = 0; i < n; ++i) {
    M.col_indices.push_back(i);
    M.values.push_back(value);
    M.row_ptr[i + 1] = i + 1;
  }
  return M;
}

static CRSMatrix BuildScatteredMatrix(int n, const Complex &value, int step = 7) {
  CRSMatrix M(n, n);
  for (int i = 0; i < n; ++i) {
    int col = (i * step) % n;
    M.col_indices.push_back(col);
    M.values.push_back(value);
    M.row_ptr[i + 1] = i + 1;
  }

  return M;
}

static CRSMatrix BuildBlockDiagonalMatrix(int n, int block_size, const Complex &value) {
  CRSMatrix M(n, n);
  for (int i = 0; i < n; ++i) {
    int block_start = (i / block_size) * block_size;
    int block_end = std::min(block_start + block_size, n);
    for (int col = block_start; col < block_end; ++col) {
      M.col_indices.push_back(col);
      M.values.push_back(value);
    }
    M.row_ptr[i + 1] = static_cast<int>(M.col_indices.size());
  }
  return M;
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, MeasureTime) {
  const int N = 1000;

  CRSMatrix A = BuildBandMatrix(N, Complex(1.0, 1.0));
  CRSMatrix B = BuildBandMatrix(N, Complex(1.0, 0.0));

  ASSERT_TRUE(A.IsValid()) << "Matrix A failed CRS validation";
  ASSERT_TRUE(B.IsValid()) << "Matrix B failed CRS validation";

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());

  auto start = std::chrono::high_resolution_clock::now();
  ASSERT_TRUE(task.Run());
  auto end = std::chrono::high_resolution_clock::now();

  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ASSERT_TRUE(task.PostProcessing());

  EXPECT_LT(duration_ms.count(), 5000) << "Multiplication took " << duration_ms.count()
                                       << " ms, which exceeds the 5 s limit";

  const CRSMatrix &C = task.GetOutput();
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, N);
  EXPECT_EQ(C.cols, N);
  EXPECT_GT(C.NonZeros(), 0u);
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, DiagonalMatrix5000x5000) {
  const int N = 5000;

  CRSMatrix A = BuildDiagonalMatrix(N, Complex(2.0, 1.0));
  CRSMatrix B = BuildDiagonalMatrix(N, Complex(1.0, -1.0));

  ASSERT_TRUE(A.IsValid());
  ASSERT_TRUE(B.IsValid());

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());

  auto start = std::chrono::high_resolution_clock::now();
  ASSERT_TRUE(task.Run());
  auto end = std::chrono::high_resolution_clock::now();

  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ASSERT_TRUE(task.PostProcessing());

  EXPECT_LT(duration_ms.count(), 500) << "Diagonal 5000x5000 took " << duration_ms.count()
                                      << " ms, expected under 500 ms";

  const CRSMatrix &C = task.GetOutput();
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, N);
  EXPECT_EQ(C.cols, N);

  EXPECT_EQ(C.NonZeros(), static_cast<size_t>(N));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, WideBandMatrix2000x2000) {
  const int N = 2000;

  CRSMatrix A = BuildBandMatrix(N, Complex(1.0, 0.5), 20);
  CRSMatrix B = BuildBandMatrix(N, Complex(0.5, 1.0), 20);

  ASSERT_TRUE(A.IsValid());
  ASSERT_TRUE(B.IsValid());

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());

  auto start = std::chrono::high_resolution_clock::now();
  ASSERT_TRUE(task.Run());
  auto end = std::chrono::high_resolution_clock::now();

  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ASSERT_TRUE(task.PostProcessing());

  EXPECT_LT(duration_ms.count(), 5000) << "Wide-band 2000x2000 took " << duration_ms.count()
                                       << " ms, exceeds 5 s limit";

  const CRSMatrix &C = task.GetOutput();
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, N);
  EXPECT_EQ(C.cols, N);
  EXPECT_GT(C.NonZeros(), 0u);
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, ScatteredMatrix3000x3000) {
  const int N = 3000;

  CRSMatrix A = BuildScatteredMatrix(N, Complex(1.0, 1.0), 7);
  CRSMatrix B = BuildScatteredMatrix(N, Complex(1.0, -1.0), 11);

  ASSERT_TRUE(A.IsValid());
  ASSERT_TRUE(B.IsValid());

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());

  auto start = std::chrono::high_resolution_clock::now();
  ASSERT_TRUE(task.Run());
  auto end = std::chrono::high_resolution_clock::now();

  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ASSERT_TRUE(task.PostProcessing());

  EXPECT_LT(duration_ms.count(), 1000) << "Scattered 3000x3000 took " << duration_ms.count()
                                       << " ms, expected under 1 s";

  const CRSMatrix &C = task.GetOutput();
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, N);
  EXPECT_EQ(C.cols, N);
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, BlockDiagonalMatrix1000x1000) {
  const int N = 1000;
  const int BLOCK_SIZE = 10;

  CRSMatrix A = BuildBlockDiagonalMatrix(N, BLOCK_SIZE, Complex(1.0, 1.0));
  CRSMatrix B = BuildBlockDiagonalMatrix(N, BLOCK_SIZE, Complex(2.0, -1.0));

  ASSERT_TRUE(A.IsValid());
  ASSERT_TRUE(B.IsValid());

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());

  auto start = std::chrono::high_resolution_clock::now();
  ASSERT_TRUE(task.Run());
  auto end = std::chrono::high_resolution_clock::now();

  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ASSERT_TRUE(task.PostProcessing());

  EXPECT_LT(duration_ms.count(), 2000) << "Block-diagonal 1000x1000 took " << duration_ms.count()
                                       << " ms, exceeds 2 s limit";

  const CRSMatrix &C = task.GetOutput();
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, N);
  EXPECT_EQ(C.cols, N);

  EXPECT_EQ(C.NonZeros(), static_cast<size_t>(N * BLOCK_SIZE));
}

TEST(VidermanASparseMatrixMultCRSComplexPerfTest, RectangularBandMatrix500x2000x500) {
  const int M = 500;
  const int K = 2000;
  const int N_out = 500;

  CRSMatrix A(M, K);
  for (int i = 0; i < M; ++i) {
    for (int offset = 0; offset < 5; ++offset) {
      int col = i * (K / M) + offset;
      if (col < K) {
        A.col_indices.push_back(col);
        A.values.push_back(Complex(1.0, 0.5));
      }
    }
    A.row_ptr[i + 1] = static_cast<int>(A.col_indices.size());
  }

  CRSMatrix B(K, N_out);
  for (int i = 0; i < K; ++i) {
    int col = i * N_out / K;
    if (col < N_out) {
      B.col_indices.push_back(col);
      B.values.push_back(Complex(0.5, 1.0));
    }
    B.row_ptr[i + 1] = static_cast<int>(B.col_indices.size());
  }

  ASSERT_TRUE(A.IsValid());
  ASSERT_TRUE(B.IsValid());

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());

  auto start = std::chrono::high_resolution_clock::now();
  ASSERT_TRUE(task.Run());
  auto end = std::chrono::high_resolution_clock::now();

  auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ASSERT_TRUE(task.PostProcessing());

  EXPECT_LT(duration_ms.count(), 2000) << "Rectangular 500x2000x500 took " << duration_ms.count()
                                       << " ms, exceeds 2 s limit";

  const CRSMatrix &C = task.GetOutput();
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, M);
  EXPECT_EQ(C.cols, N_out);
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
