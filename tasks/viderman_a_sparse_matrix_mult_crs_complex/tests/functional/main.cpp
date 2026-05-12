#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {
namespace {
constexpr double kTestTol = 1e-12;

bool ComplexNear(const Complex &lhs, const Complex &rhs, double tol = kTestTol) {
  return std::abs(lhs.real() - rhs.real()) <= tol && std::abs(lhs.imag() - rhs.imag()) <= tol;
}

bool CrsEqual(const CRSMatrix &expected, const CRSMatrix &actual, double tol = kTestTol) {
  if (expected.rows != actual.rows || expected.cols != actual.cols) {
    return false;
  }
  if (expected.row_ptr != actual.row_ptr || expected.col_indices != actual.col_indices) {
    return false;
  }
  if (expected.values.size() != actual.values.size()) {
    return false;
  }
  for (std::size_t i = 0; i < expected.values.size(); ++i) {
    if (!ComplexNear(expected.values[i], actual.values[i], tol)) {
      return false;
    }
  }
  return true;
}

void CompareCrsMatrices(const CRSMatrix &expected, const CRSMatrix &actual) {
  EXPECT_TRUE(CrsEqual(expected, actual));
}

std::vector<std::vector<Complex>> ToDense(const CRSMatrix &m) {
  std::vector<std::vector<Complex>> dense(m.rows, std::vector<Complex>(m.cols, {0.0, 0.0}));
  for (int i = 0; i < m.rows; ++i) {
    for (int j = m.row_ptr[i]; j < m.row_ptr[i + 1]; ++j) {
      dense[i][m.col_indices[j]] = m.values[j];
    }
  }
  return dense;
}

bool DenseEqual(const std::vector<std::vector<Complex>> &expected, const CRSMatrix &actual, double tol = kTestTol) {
  const int rows = static_cast<int>(expected.size());
  const int cols = rows > 0 ? static_cast<int>(expected[0].size()) : 0;
  if (actual.rows != rows || actual.cols != cols) {
    return false;
  }
  const auto dense = ToDense(actual);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (!ComplexNear(dense[i][j], expected[i][j], tol)) {
        return false;
      }
    }
  }
  return true;
}

void CompareDense(const std::vector<std::vector<Complex>> &expected, const CRSMatrix &actual, double tol = kTestTol) {
  EXPECT_TRUE(DenseEqual(expected, actual, tol));
}

bool Dense2x2Equal(const CRSMatrix &lhs, const CRSMatrix &rhs, double tol = 1e-11) {
  const auto l = ToDense(lhs);
  const auto r = ToDense(rhs);
  if (l.size() != 2 || r.size() != 2 || l[0].size() != 2 || r[0].size() != 2) {
    return false;
  }
  return ComplexNear(l[0][0], r[0][0], tol) && ComplexNear(l[0][1], r[0][1], tol) &&
         ComplexNear(l[1][0], r[1][0], tol) && ComplexNear(l[1][1], r[1][1], tol);
}

void Compare2x2Dense(const CRSMatrix &lhs, const CRSMatrix &rhs, double tol = 1e-11) {
  EXPECT_TRUE(Dense2x2Equal(lhs, rhs, tol));
}

bool CheckPartialCancellation(const CRSMatrix &c) {
  return c.rows == 2 && c.cols == 1 && c.row_ptr == std::vector<int>({0, 0, 1}) && c.values.size() == 1 &&
         ComplexNear(c.values[0], Complex(0.0, 2.0));
}

bool CheckCornerElements5x5(const CRSMatrix &c) {
  const auto dense = ToDense(c);
  return c.rows == 5 && c.cols == 5 && c.NonZeros() == 2U && ComplexNear(dense[0][0], Complex(2.0, 0.0)) &&
         ComplexNear(dense[4][4], Complex(1.0, 0.0));
}

bool CheckDenseRowTimesIdentity(const CRSMatrix &c) {
  const auto dense = ToDense(c);
  return c.NonZeros() == 4U && ComplexNear(dense[0][0], Complex(1.0, 1.0)) &&
         ComplexNear(dense[0][1], Complex(2.0, 0.0)) && ComplexNear(dense[0][2], Complex(3.0, -1.0)) &&
         ComplexNear(dense[0][3], Complex(0.0, 4.0));
}

bool IsShapeAndEmpty(const CRSMatrix &c, int rows, int cols) {
  return c.rows == rows && c.cols == cols && c.values.empty();
}

bool IsShapeAndValid(const CRSMatrix &c, int rows, int cols) {
  return c.rows == rows && c.cols == cols && c.IsValid();
}

bool IsSortedInRows(const CRSMatrix &c) {
  for (int i = 0; i < c.rows; ++i) {
    for (int j = c.row_ptr[i]; j < c.row_ptr[i + 1] - 1; ++j) {
      if (c.col_indices[j] >= c.col_indices[j + 1]) {
        return false;
      }
    }
  }
  return true;
}

bool HasNoExplicitZeros(const CRSMatrix &c) {
  return std::ranges::all_of(c.values, [](const Complex &v) { return std::abs(v) > kEpsilon; });
}

CRSMatrix RunTask(const CRSMatrix &a, const CRSMatrix &b, bool expect_valid = true) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(a, b));
  if (!expect_valid) {
    EXPECT_FALSE(task.Validation());
    return CRSMatrix{};
  }
  const bool ok = task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
  EXPECT_TRUE(ok);
  return ok ? task.GetOutput() : CRSMatrix{};
}
}  // namespace

TEST(VidermanValidation, IncompatibleDimensions) {
  RunTask(CRSMatrix(2, 3), CRSMatrix(4, 5), false);
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

  CompareCrsMatrices(expected, RunTask(a, b));
}

TEST(VidermanEdgeCases, BothZeroMatrices) {
  EXPECT_TRUE(IsShapeAndValid(RunTask(CRSMatrix(3, 4), CRSMatrix(4, 5)), 3, 5));
}

TEST(VidermanEdgeCases, ZeroANonzeroB) {
  CRSMatrix a(2, 3);
  CRSMatrix b(3, 2);
  b.row_ptr = {0, 1, 2, 3};
  b.col_indices = {0, 1, 0};
  b.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};
  EXPECT_TRUE(IsShapeAndEmpty(RunTask(a, b), 2, 2));
}

TEST(VidermanEdgeCases, NonzeroAZeroB) {
  CRSMatrix a(2, 3);
  a.row_ptr = {0, 2, 3};
  a.col_indices = {0, 2, 1};
  a.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};
  EXPECT_TRUE(IsShapeAndEmpty(RunTask(a, CRSMatrix(3, 4)), 2, 4));
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
  CompareCrsMatrices(expected, RunTask(a, b));
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

  const std::vector<std::vector<Complex>> expected = {{Complex(3, 3), Complex(-1, 1)}, {Complex(6, 0), Complex(0, 2)}};
  CompareDense(expected, RunTask(a, b));
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
  const std::vector<Complex> a_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  const std::vector<Complex> b_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      expected[i][j] = a_vals[i] * b_vals[j];
    }
  }
  CompareDense(expected, RunTask(a, b));
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
  CompareCrsMatrices(expected, RunTask(a, b));
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
  CompareCrsMatrices(expected, RunTask(a, b));
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
  CompareCrsMatrices(expected, RunTask(a, b));
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
  EXPECT_TRUE(IsShapeAndEmpty(RunTask(a, b), 1, 1));
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
  EXPECT_TRUE(CheckPartialCancellation(RunTask(a, b)));
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
  CompareCrsMatrices(a, RunTask(a, i));
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
  CompareCrsMatrices(a, RunTask(i, a));
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
  Compare2x2Dense(RunTask(RunTask(a, b), c), RunTask(a, RunTask(b, c)));
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
  CompareCrsMatrices(minus_i, RunTask(i_i, i_i));
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
  CompareCrsMatrices(i, RunTask(p, pt));
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
  EXPECT_TRUE(IsShapeAndValid(RunTask(a, b), 4, 5));
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
  EXPECT_TRUE(IsSortedInRows(RunTask(a, b)));
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
  const CRSMatrix c = RunTask(a, b);
  EXPECT_TRUE(!c.row_ptr.empty() && c.row_ptr[0] == 0);
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
  const CRSMatrix c = RunTask(a, b);
  EXPECT_TRUE(static_cast<std::size_t>(c.row_ptr[c.rows]) == c.values.size());
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
  EXPECT_TRUE(HasNoExplicitZeros(RunTask(a, b)));
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

  const std::vector<std::vector<Complex>> expected = {{Complex(6, 3), Complex(1, 1), Complex(0, 0), Complex(0, 0)},
                                                      {Complex(0, 0), Complex(0, 0), Complex(6, 0), Complex(3, 3)}};
  CompareDense(expected, RunTask(a, b));
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
  CompareCrsMatrices(expected, RunTask(a, b));
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
  EXPECT_TRUE(CheckCornerElements5x5(RunTask(a, b)));
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
  EXPECT_TRUE(CheckDenseRowTimesIdentity(RunTask(a, i)));
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
  CompareCrsMatrices(minus_i, RunTask(j, j));
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
