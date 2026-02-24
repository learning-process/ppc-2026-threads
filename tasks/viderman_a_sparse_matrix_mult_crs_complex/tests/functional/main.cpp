#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {
constexpr double kTestTol = 1e-12;

bool complex_near(const Complex &lhs, const Complex &rhs, double tol = kTestTol) {
  return std::abs(lhs.real() - rhs.real()) <= tol && std::abs(lhs.imag() - rhs.imag()) <= tol;
}

bool crs_equal(const CRSMatrix &expected, const CRSMatrix &actual, double tol = kTestTol) {
  if (expected.rows != actual.rows || expected.cols != actual.cols) {
    return false;
  }
  if (expected.row_ptr != actual.row_ptr || expected.col_indices != actual.col_indices) {
    return false;
  }
  if (expected.values.size() != actual.values.size()) {
    return false;
  }
  for (size_t i = 0; i < expected.values.size(); ++i) {
    if (!complex_near(expected.values[i], actual.values[i], tol)) {
      return false;
    }
  }
  return true;
}

void compare_crs_matrices(const CRSMatrix &expected, const CRSMatrix &actual) {
  EXPECT_TRUE(crs_equal(expected, actual));
}

std::vector<std::vector<Complex>> to_dense(const CRSMatrix &m) {
  std::vector<std::vector<Complex>> d(m.rows, std::vector<Complex>(m.cols, {0.0, 0.0}));
  for (int i = 0; i < m.rows; ++i) {
    for (int j = m.row_ptr[i]; j < m.row_ptr[i + 1]; ++j) {
      d[i][m.col_indices[j]] = m.values[j];
    }
  }
  return d;
}

void compare_dense(const std::vector<std::vector<Complex>> &expected, const CRSMatrix &actual, double tol = 1e-12) {
  const int rows = static_cast<int>(expected.size());
  const int cols = rows > 0 ? static_cast<int>(expected[0].size()) : 0;
  bool ok = actual.rows == rows && actual.cols == cols;
  if (ok) {
    const auto d = to_dense(actual);
    const int total = rows * cols;
    for (int idx = 0; idx < total; ++idx) {
      const int i = idx / cols;
      const int j = idx % cols;
      if (!complex_near(d[i][j], expected[i][j], tol)) {
        ok = false;
        break;
      }
    }
  }
  EXPECT_TRUE(ok);
}

void compare_2x2_dense(const CRSMatrix &lhs, const CRSMatrix &rhs, double tol = 1e-11) {
  const auto l = to_dense(lhs);
  const auto r = to_dense(rhs);
  bool ok = l.size() == r.size() && l.size() == 2;
  for (int i = 0; ok && i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      if (!complex_near(l[i][j], r[i][j], tol)) {
        ok = false;
        break;
      }
    }
  }
  EXPECT_TRUE(ok);
}

void check_partial_cancellation(const CRSMatrix &c) {
  const bool ok = c.rows == 2 && c.cols == 1 && c.row_ptr.size() == 3 && c.row_ptr[0] == 0 && c.row_ptr[1] == 0 &&
                  c.row_ptr[2] == 1 && c.values.size() == 1 && complex_near(c.values[0], Complex(0.0, 2.0));
  EXPECT_TRUE(ok);
}

void check_corner_elements_5x5(const CRSMatrix &c) {
  const auto d = to_dense(c);
  const bool ok = c.rows == 5 && c.cols == 5 && c.non_zeros() == 2U && complex_near(d[0][0], Complex(2.0, 0.0)) &&
                  complex_near(d[4][4], Complex(1.0, 0.0));
  EXPECT_TRUE(ok);
}

void check_dense_row_times_identity(const CRSMatrix &c) {
  const auto d = to_dense(c);
  const bool ok = c.non_zeros() == 4U && complex_near(d[0][0], Complex(1.0, 1.0)) &&
                  complex_near(d[0][1], Complex(2.0, 0.0)) && complex_near(d[0][2], Complex(3.0, -1.0)) &&
                  complex_near(d[0][3], Complex(0.0, 4.0));
  EXPECT_TRUE(ok);
}

CRSMatrix run_task(const CRSMatrix &a, const CRSMatrix &b, bool expect_valid = true) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  if (!expect_valid) {
    EXPECT_FALSE(task.Validation());
    return CRSMatrix{};
  }
  const bool ok = task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
  EXPECT_TRUE(ok);
  if (!ok) {
    return CRSMatrix{};
  }
  return task.GetOutput();
}
}  // namespace

TEST(VidermanValidation, IncompatibleDimensions) {
  run_task(CRSMatrix(2, 3), CRSMatrix(4, 5), false);
}

TEST(VidermanValidation, NegativeColIndex) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, -1};
  a.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, ColIndexOutOfRange) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 5};
  a.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix b(2, 3);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, UnsortedColIndices) {
  CRSMatrix a(1, 3);
  a.row_ptr = {0, 3};
  a.col_indices = {2, 0, 1};
  a.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, NonMonotonicRowPtr) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 1};
  a.col_indices = {0, 1, 0};
  a.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, WrongRowPtrSize) {
  CRSMatrix a;
  a.rows = 3;
  a.cols = 3;
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, ColIndicesValuesSizeMismatch) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanEdgeCases, SingleElement) {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3.0, 4.0)};

  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(1.0, -2.0)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(11.0, -2.0)};

  compare_crs_matrices(expected, run_task(a, b));
}

TEST(VidermanEdgeCases, BothZeroMatrices) {
  CRSMatrix a(3, 4);
  CRSMatrix b(4, 5);

  CRSMatrix c = run_task(a, b);
  EXPECT_EQ(c.rows, 3);
  EXPECT_EQ(c.cols, 5);
  EXPECT_TRUE(c.values.empty());
  EXPECT_TRUE(c.is_valid());
}

TEST(VidermanEdgeCases, ZeroANonzeroB) {
  CRSMatrix a(2, 3);
  CRSMatrix b(3, 2);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 0};
  b.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix c = run_task(a, b);
  EXPECT_EQ(c.rows, 2);
  EXPECT_EQ(c.cols, 2);
  EXPECT_TRUE(c.values.empty());
}

TEST(VidermanEdgeCases, NonzeroAZeroB) {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};

  CRSMatrix b(3, 4);

  CRSMatrix c = run_task(a, b);
  EXPECT_EQ(c.rows, 2);
  EXPECT_EQ(c.cols, 4);
  EXPECT_TRUE(c.values.empty());
}

TEST(VidermanEdgeCases, RowVectorTimesColVector) {
  CRSMatrix a(1, 3);
  a.row_ptr = {0, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};

  CRSMatrix b(3, 1);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 0, 0};
  b.values = {Complex(1, 0), Complex(0, 1), Complex(1, 1)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(5.0, 5.0)};

  compare_crs_matrices(expected, run_task(a, b));
}

TEST(VidermanEdgeCases, ColVectorTimesRowVector) {
  CRSMatrix a(2, 1);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 0};
  a.values = {Complex(1, 1), Complex(2, 0)};

  CRSMatrix b(1, 2);
  b.row_ptr = {0, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(3, 0), Complex(0, 1)};

  std::vector<std::vector<Complex>> expected = {{Complex(3, 3), Complex(-1, 1)}, {Complex(6, 0), Complex(0, 2)}};
  compare_dense(expected, run_task(a, b));
}

TEST(VidermanEdgeCases, TallSkinnyMatrix) {
  CRSMatrix a(4, 1);
  a.row_ptr = {0, 1, 2, 3, 4};
  a.col_indices = {0, 0, 0, 0};
  a.values = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};

  CRSMatrix b(1, 4);
  b.row_ptr = {0, 4};
  b.col_indices = {0, 1, 2, 3};
  b.values = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};

  std::vector<std::vector<Complex>> expected(4, std::vector<Complex>(4));
  std::vector<Complex> a_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  std::vector<Complex> b_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      expected[i][j] = a_vals[i] * b_vals[j];
    }
  }

  compare_dense(expected, run_task(a, b));
}

TEST(VidermanComplexArithmetic, DiagonalMultiplication) {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 1, 2, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 1), Complex(2, 2), Complex(3, 3)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 2};
  b.values = {Complex(4, 0), Complex(5, 0), Complex(6, 0)};

  CRSMatrix expected(3, 3);
  expected.row_ptr = {0, 1, 2, 3};
  expected.col_indices = {0, 1, 2};
  expected.values = {Complex(4, 4), Complex(10, 10), Complex(18, 18)};

  compare_crs_matrices(expected, run_task(a, b));
}

TEST(VidermanComplexArithmetic, PureImaginarySquared) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix expected(2, 2);
  expected.row_ptr = {0, 1, 2};
  expected.col_indices = {0, 1};
  expected.values = {Complex(-1, 0), Complex(-1, 0)};

  compare_crs_matrices(expected, run_task(a, b));
}

TEST(VidermanComplexArithmetic, ConjugateProduct) {
  CRSMatrix a(1, 1);
  a.row_ptr = {0, 1};
  a.col_indices = {0};
  a.values = {Complex(3, 4)};

  CRSMatrix b(1, 1);
  b.row_ptr = {0, 1};
  b.col_indices = {0};
  b.values = {Complex(3, -4)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(25, 0)};

  compare_crs_matrices(expected, run_task(a, b));
}

TEST(VidermanComplexArithmetic, CancellationToZero) {
  CRSMatrix a(1, 2);
  a.row_ptr = {0, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(-1, 0)};

  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix c = run_task(a, b);
  EXPECT_EQ(c.rows, 1);
  EXPECT_EQ(c.cols, 1);
  EXPECT_TRUE(c.values.empty()) << "Ожидалась нулевая матрица, но нашлось " << c.values.size() << " элемент(ов)";
}

TEST(VidermanComplexArithmetic, PartialCancellation) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 4};
  a.col_indices = {0, 1, 0, 1};
  a.values = {Complex(1, 0), Complex(-1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(0, 1), Complex(0, 1)};

  const CRSMatrix c = run_task(a, b);
  check_partial_cancellation(c);
}

TEST(VidermanAlgebraic, RightIdentity) {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 2, 3, 4};
  a.col_indices = {0, 2, 1, 0};
  a.values = {Complex(1, 2), Complex(3, 0), Complex(0, 1), Complex(5, -1)};

  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  compare_crs_matrices(a, run_task(a, i));
}

TEST(VidermanAlgebraic, LeftIdentity) {
  CRSMatrix a(3, 3);
  a.row_ptr = {0, 2, 3, 4};
  a.col_indices = {0, 2, 1, 0};
  a.values = {Complex(1, 2), Complex(3, 0), Complex(0, 1), Complex(5, -1)};

  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  compare_crs_matrices(a, run_task(i, a));
}

TEST(VidermanAlgebraic, Associativity) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 1, 0};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(0, 1)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(3, 0), Complex(0, 2)};

  CRSMatrix c(2, 2);
  c.row_ptr = {0, 2, 2};
  c.col_indices = {0, 1};
  c.values = {Complex(1, 0), Complex(1, 1)};

  CRSMatrix ab = run_task(a, b);
  CRSMatrix abc_left = run_task(ab, c);

  CRSMatrix bc = run_task(b, c);
  CRSMatrix abc_right = run_task(a, bc);

  compare_2x2_dense(abc_left, abc_right);
}

TEST(VidermanAlgebraic, SquareOfScaledIdentity) {
  CRSMatrix i_i(2, 2);
  i_i.row_ptr = {0, 1, 2};
  i_i.col_indices = {0, 1};
  i_i.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix minus_i(2, 2);
  minus_i.row_ptr = {0, 1, 2};
  minus_i.col_indices = {0, 1};
  minus_i.values = {Complex(-1, 0), Complex(-1, 0)};

  compare_crs_matrices(minus_i, run_task(i_i, i_i));
}

TEST(VidermanAlgebraic, PermutationTimesTranspose) {
  CRSMatrix p(3, 3);
  p.row_ptr = {0, 1, 2, 3};
  p.col_indices = {1, 2, 0};
  p.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix pt(3, 3);
  pt.row_ptr = {0, 1, 2, 3};
  pt.col_indices = {2, 0, 1};
  pt.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix i(3, 3);
  i.row_ptr = {0, 1, 2, 3};
  i.col_indices = {0, 1, 2};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  compare_crs_matrices(i, run_task(p, pt));
}

TEST(VidermanStructural, OutputIsValidCRS) {
  CRSMatrix a(4, 3);
  a.row_ptr = {0, 2, 3, 3, 4};
  a.col_indices = {0, 2, 1, 2};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(0, 3), Complex(-1, 1)};

  CRSMatrix b(3, 5);
  b.row_ptr = {0, 2, 3, 5};
  b.col_indices = {0, 3, 2, 1, 4};
  b.values = {Complex(1, 0), Complex(2, 1), Complex(0, 1), Complex(3, 0), Complex(1, -1)};

  CRSMatrix c = run_task(a, b);
  EXPECT_TRUE(c.is_valid());
  EXPECT_EQ(c.rows, 4);
  EXPECT_EQ(c.cols, 5);
}

TEST(VidermanStructural, ColIndicesSortedInOutput) {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 3, 3};
  a.col_indices = {0, 1, 2};
  a.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(3, 3);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {2, 0, 1};
  b.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix c = run_task(a, b);
  for (int i = 0; i < c.rows; ++i) {
    for (int j = c.row_ptr[i]; j < c.row_ptr[i + 1] - 1; ++j) {
      EXPECT_LT(c.col_indices[j], c.col_indices[j + 1])
          << "Строка " << i << ": col_indices не отсортированы на позиции " << j;
    }
  }
}

TEST(VidermanStructural, RowPtrStartsAtZero) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 1};
  b.values = {Complex(3, 0), Complex(4, 0)};

  CRSMatrix c = run_task(a, b);
  ASSERT_FALSE(c.row_ptr.empty());
  EXPECT_EQ(c.row_ptr[0], 0);
}

TEST(VidermanStructural, RowPtrLastEqualsNNZ) {
  CRSMatrix a(3, 2);
  a.row_ptr = {0, 2, 2, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 1), Complex(2, 0)};

  CRSMatrix b(2, 3);
  b.row_ptr = {0, 2, 3};
  b.col_indices = {0, 2, 1};
  b.values = {Complex(1, 0), Complex(0, 1), Complex(3, 0)};

  CRSMatrix c = run_task(a, b);
  EXPECT_EQ(static_cast<size_t>(c.row_ptr[c.rows]), c.values.size());
}

TEST(VidermanStructural, NoExplicitZerosInOutput) {
  CRSMatrix a(2, 2);
  a.row_ptr = {0, 1, 2};
  a.col_indices = {0, 1};
  a.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix b(2, 2);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {1, 0};
  b.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix c = run_task(a, b);
  for (const auto &v : c.values) {
    EXPECT_GT(std::abs(v), 1e-14) << "В values присутствует явный ноль";
  }
}

TEST(VidermanNonTrivial, RectangularWithAccumulation) {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {Complex(1, 0), Complex(2, 1), Complex(3, 0)};

  CRSMatrix b(3, 4);
  b.row_ptr = {0, 1, 3, 4};
  b.col_indices = {1, 2, 3, 0};
  b.values = {Complex(1, 1), Complex(2, 0), Complex(1, 1), Complex(3, 0)};

  std::vector<std::vector<Complex>> expected = {{Complex(6, 3), Complex(1, 1), Complex(0, 0), Complex(0, 0)},
                                                {Complex(0, 0), Complex(0, 0), Complex(6, 0), Complex(3, 3)}};
  compare_dense(expected, run_task(a, b));
}

TEST(VidermanNonTrivial, MultipleRowsContributeToSameColumn) {
  CRSMatrix a(3, 2);
  a.row_ptr = {0, 2, 4, 6};
  a.col_indices = {0, 1, 0, 1, 0, 1};
  a.values = {Complex(1, 0), Complex(0, 1), Complex(2, 0), Complex(0, 2), Complex(3, 0), Complex(0, 3)};

  CRSMatrix b(2, 1);
  b.row_ptr = {0, 1, 2};
  b.col_indices = {0, 0};
  b.values = {Complex(1, 1), Complex(1, -1)};

  CRSMatrix expected(3, 1);
  expected.row_ptr = {0, 1, 2, 3};
  expected.col_indices = {0, 0, 0};
  expected.values = {Complex(2, 2), Complex(4, 4), Complex(6, 6)};

  compare_crs_matrices(expected, run_task(a, b));
}

TEST(VidermanNonTrivial, CornerElementsOnly5x5) {
  CRSMatrix a(5, 5);
  a.row_ptr = {0, 1, 1, 1, 1, 2};
  a.col_indices = {0, 4};
  a.values = {Complex(1, 0), Complex(0, 1)};

  CRSMatrix b(5, 5);
  b.row_ptr = {0, 1, 1, 1, 1, 2};
  b.col_indices = {0, 4};
  b.values = {Complex(2, 0), Complex(0, -1)};

  const CRSMatrix c = run_task(a, b);
  check_corner_elements_5x5(c);
}

TEST(VidermanNonTrivial, DenseRowTimesIdentity) {
  CRSMatrix a(1, 4);
  a.row_ptr = {0, 4};
  a.col_indices = {0, 1, 2, 3};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1), Complex(0, 4)};

  CRSMatrix i(4, 4);
  i.row_ptr = {0, 1, 2, 3, 4};
  i.col_indices = {0, 1, 2, 3};
  i.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  const CRSMatrix c = run_task(a, i);
  check_dense_row_times_identity(c);
}

TEST(VidermanNonTrivial, MatrixSquaredKnownResult) {
  CRSMatrix j(2, 2);
  j.row_ptr = {0, 1, 2};
  j.col_indices = {1, 0};
  j.values = {Complex(1, 0), Complex(-1, 0)};

  CRSMatrix minus_i(2, 2);
  minus_i.row_ptr = {0, 1, 2};
  minus_i.col_indices = {0, 1};
  minus_i.values = {Complex(-1, 0), Complex(-1, 0)};

  compare_crs_matrices(minus_i, run_task(j, j));
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
