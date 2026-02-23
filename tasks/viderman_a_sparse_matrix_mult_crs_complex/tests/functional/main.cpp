#include <gtest/gtest.h>

#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/seq/include/ops_seq.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

void CompareCRSMatrices(const CRSMatrix &expected, const CRSMatrix &actual) {
  EXPECT_EQ(expected.rows, actual.rows);
  EXPECT_EQ(expected.cols, actual.cols);
  EXPECT_EQ(expected.row_ptr, actual.row_ptr);
  EXPECT_EQ(expected.col_indices, actual.col_indices);
  ASSERT_EQ(expected.values.size(), actual.values.size());
  for (size_t i = 0; i < expected.values.size(); ++i) {
    EXPECT_NEAR(expected.values[i].real(), actual.values[i].real(), 1e-12);
    EXPECT_NEAR(expected.values[i].imag(), actual.values[i].imag(), 1e-12);
  }
}

std::vector<std::vector<Complex>> ToDense(const CRSMatrix &M) {
  std::vector<std::vector<Complex>> D(M.rows, std::vector<Complex>(M.cols, {0.0, 0.0}));
  for (int i = 0; i < M.rows; ++i) {
    for (int j = M.row_ptr[i]; j < M.row_ptr[i + 1]; ++j) {
      D[i][M.col_indices[j]] = M.values[j];
    }
  }
  return D;
}

void CompareDense(const std::vector<std::vector<Complex>> &expected, const CRSMatrix &actual, double tol = 1e-12) {
  int rows = static_cast<int>(expected.size());
  int cols = rows > 0 ? static_cast<int>(expected[0].size()) : 0;
  ASSERT_EQ(actual.rows, rows);
  ASSERT_EQ(actual.cols, cols);
  auto D = ToDense(actual);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      EXPECT_NEAR(D[i][j].real(), expected[i][j].real(), tol) << "  at [" << i << "][" << j << "] real part";
      EXPECT_NEAR(D[i][j].imag(), expected[i][j].imag(), tol) << "  at [" << i << "][" << j << "] imag part";
    }
  }
}

CRSMatrix RunTask(const CRSMatrix &A, const CRSMatrix &B, bool expect_valid = true) {
  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  if (!expect_valid) {
    EXPECT_FALSE(task.Validation());
    return CRSMatrix{};
  }
  EXPECT_TRUE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());
  EXPECT_TRUE(task.Run());
  EXPECT_TRUE(task.PostProcessing());
  return task.GetOutput();
}

TEST(VidermanValidation, IncompatibleDimensions) {
  RunTask(CRSMatrix(2, 3), CRSMatrix(4, 5), false);
}

TEST(VidermanValidation, NegativeColIndex) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, -1};
  A.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, ColIndexOutOfRange) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 5};
  A.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix B(2, 3);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, UnsortedColIndices) {
  CRSMatrix A(1, 3);
  A.row_ptr = {0, 3};
  A.col_indices = {2, 0, 1};
  A.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix B(3, 3);
  B.row_ptr = {0, 1, 2, 3};
  B.col_indices = {0, 1, 2};
  B.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, NonMonotonicRowPtr) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 2, 1};
  A.col_indices = {0, 1, 0};
  A.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, WrongRowPtrSize) {
  CRSMatrix A;
  A.rows = 3;
  A.cols = 3;
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix B(3, 3);
  B.row_ptr = {0, 1, 2, 3};
  B.col_indices = {0, 1, 2};
  B.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanValidation, ColIndicesValuesSizeMismatch) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(1, 0)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(1, 0), Complex(1, 0)};

  VidermanASparseMatrixMultCRSComplexSEQ task(std::make_tuple(A, B));
  EXPECT_FALSE(task.Validation());
}

TEST(VidermanEdgeCases, SingleElement) {
  CRSMatrix A(1, 1);
  A.row_ptr = {0, 1};
  A.col_indices = {0};
  A.values = {Complex(3.0, 4.0)};

  CRSMatrix B(1, 1);
  B.row_ptr = {0, 1};
  B.col_indices = {0};
  B.values = {Complex(1.0, -2.0)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(11.0, -2.0)};

  CompareCRSMatrices(expected, RunTask(A, B));
}

TEST(VidermanEdgeCases, BothZeroMatrices) {
  CRSMatrix A(3, 4);
  CRSMatrix B(4, 5);

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(C.rows, 3);
  EXPECT_EQ(C.cols, 5);
  EXPECT_TRUE(C.values.empty());
  EXPECT_TRUE(C.IsValid());
}

TEST(VidermanEdgeCases, ZeroANonzeroB) {
  CRSMatrix A(2, 3);
  CRSMatrix B(3, 2);
  B.row_ptr = {0, 1, 2, 3};
  B.col_indices = {0, 1, 0};
  B.values = {Complex(1, 0), Complex(2, 0), Complex(3, 0)};

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(C.rows, 2);
  EXPECT_EQ(C.cols, 2);
  EXPECT_TRUE(C.values.empty());
}

TEST(VidermanEdgeCases, NonzeroAZeroB) {
  CRSMatrix A(2, 3);
  A.row_ptr = {0, 2, 3};
  A.col_indices = {0, 2, 1};
  A.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};

  CRSMatrix B(3, 4);

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(C.rows, 2);
  EXPECT_EQ(C.cols, 4);
  EXPECT_TRUE(C.values.empty());
}

TEST(VidermanEdgeCases, RowVectorTimesColVector) {
  CRSMatrix A(1, 3);
  A.row_ptr = {0, 3};
  A.col_indices = {0, 1, 2};
  A.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1)};

  CRSMatrix B(3, 1);
  B.row_ptr = {0, 1, 2, 3};
  B.col_indices = {0, 0, 0};
  B.values = {Complex(1, 0), Complex(0, 1), Complex(1, 1)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(5.0, 5.0)};

  CompareCRSMatrices(expected, RunTask(A, B));
}

TEST(VidermanEdgeCases, ColVectorTimesRowVector) {
  CRSMatrix A(2, 1);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 0};
  A.values = {Complex(1, 1), Complex(2, 0)};

  CRSMatrix B(1, 2);
  B.row_ptr = {0, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(3, 0), Complex(0, 1)};

  std::vector<std::vector<Complex>> expected = {{Complex(3, 3), Complex(-1, 1)}, {Complex(6, 0), Complex(0, 2)}};
  CompareDense(expected, RunTask(A, B));
}

TEST(VidermanEdgeCases, TallSkinnyMatrix) {
  CRSMatrix A(4, 1);
  A.row_ptr = {0, 1, 2, 3, 4};
  A.col_indices = {0, 0, 0, 0};
  A.values = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};

  CRSMatrix B(1, 4);
  B.row_ptr = {0, 4};
  B.col_indices = {0, 1, 2, 3};
  B.values = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};

  std::vector<std::vector<Complex>> expected(4, std::vector<Complex>(4));
  std::vector<Complex> a_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  std::vector<Complex> b_vals = {Complex(1, 0), Complex(0, 1), Complex(-1, 0), Complex(0, -1)};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      expected[i][j] = a_vals[i] * b_vals[j];
    }
  }

  CompareDense(expected, RunTask(A, B));
}

TEST(VidermanComplexArithmetic, DiagonalMultiplication) {
  CRSMatrix A(3, 3);
  A.row_ptr = {0, 1, 2, 3};
  A.col_indices = {0, 1, 2};
  A.values = {Complex(1, 1), Complex(2, 2), Complex(3, 3)};

  CRSMatrix B(3, 3);
  B.row_ptr = {0, 1, 2, 3};
  B.col_indices = {0, 1, 2};
  B.values = {Complex(4, 0), Complex(5, 0), Complex(6, 0)};

  CRSMatrix expected(3, 3);
  expected.row_ptr = {0, 1, 2, 3};
  expected.col_indices = {0, 1, 2};
  expected.values = {Complex(4, 4), Complex(10, 10), Complex(18, 18)};

  CompareCRSMatrices(expected, RunTask(A, B));
}

TEST(VidermanComplexArithmetic, PureImaginarySquared) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix expected(2, 2);
  expected.row_ptr = {0, 1, 2};
  expected.col_indices = {0, 1};
  expected.values = {Complex(-1, 0), Complex(-1, 0)};

  CompareCRSMatrices(expected, RunTask(A, B));
}

TEST(VidermanComplexArithmetic, ConjugateProduct) {
  CRSMatrix A(1, 1);
  A.row_ptr = {0, 1};
  A.col_indices = {0};
  A.values = {Complex(3, 4)};

  CRSMatrix B(1, 1);
  B.row_ptr = {0, 1};
  B.col_indices = {0};
  B.values = {Complex(3, -4)};

  CRSMatrix expected(1, 1);
  expected.row_ptr = {0, 1};
  expected.col_indices = {0};
  expected.values = {Complex(25, 0)};

  CompareCRSMatrices(expected, RunTask(A, B));
}

TEST(VidermanComplexArithmetic, CancellationToZero) {
  CRSMatrix A(1, 2);
  A.row_ptr = {0, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(1, 0), Complex(-1, 0)};

  CRSMatrix B(2, 1);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 0};
  B.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(C.rows, 1);
  EXPECT_EQ(C.cols, 1);
  EXPECT_TRUE(C.values.empty()) << "Ожидалась нулевая матрица, но нашлось " << C.values.size() << " элемент(ов)";
}

TEST(VidermanComplexArithmetic, PartialCancellation) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 2, 4};
  A.col_indices = {0, 1, 0, 1};
  A.values = {Complex(1, 0), Complex(-1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix B(2, 1);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 0};
  B.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(C.rows, 2);
  EXPECT_EQ(C.cols, 1);

  EXPECT_EQ(C.row_ptr[0], 0);
  EXPECT_EQ(C.row_ptr[1], 0) << "Строка 0 должна быть нулевой после компенсации";

  ASSERT_EQ(C.row_ptr[2], 1) << "Строка 1 должна содержать один ненулевой элемент";
  EXPECT_NEAR(C.values[0].real(), 0.0, 1e-12);
  EXPECT_NEAR(C.values[0].imag(), 2.0, 1e-12);
}

TEST(VidermanAlgebraic, RightIdentity) {
  CRSMatrix A(3, 3);
  A.row_ptr = {0, 2, 3, 4};
  A.col_indices = {0, 2, 1, 0};
  A.values = {Complex(1, 2), Complex(3, 0), Complex(0, 1), Complex(5, -1)};

  CRSMatrix I(3, 3);
  I.row_ptr = {0, 1, 2, 3};
  I.col_indices = {0, 1, 2};
  I.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CompareCRSMatrices(A, RunTask(A, I));
}

TEST(VidermanAlgebraic, LeftIdentity) {
  CRSMatrix A(3, 3);
  A.row_ptr = {0, 2, 3, 4};
  A.col_indices = {0, 2, 1, 0};
  A.values = {Complex(1, 2), Complex(3, 0), Complex(0, 1), Complex(5, -1)};

  CRSMatrix I(3, 3);
  I.row_ptr = {0, 1, 2, 3};
  I.col_indices = {0, 1, 2};
  I.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CompareCRSMatrices(A, RunTask(I, A));
}

TEST(VidermanAlgebraic, Associativity) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 2, 3};
  A.col_indices = {0, 1, 0};
  A.values = {Complex(1, 1), Complex(2, 0), Complex(0, 1)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(3, 0), Complex(0, 2)};

  CRSMatrix C(2, 2);
  C.row_ptr = {0, 2, 2};
  C.col_indices = {0, 1};
  C.values = {Complex(1, 0), Complex(1, 1)};

  CRSMatrix AB = RunTask(A, B);
  CRSMatrix ABC_left = RunTask(AB, C);

  CRSMatrix BC = RunTask(B, C);
  CRSMatrix ABC_right = RunTask(A, BC);

  auto D_left = ToDense(ABC_left);
  auto D_right = ToDense(ABC_right);
  ASSERT_EQ(D_left.size(), D_right.size());
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_NEAR(D_left[i][j].real(), D_right[i][j].real(), 1e-11) << "real at [" << i << "][" << j << "]";
      EXPECT_NEAR(D_left[i][j].imag(), D_right[i][j].imag(), 1e-11) << "imag at [" << i << "][" << j << "]";
    }
  }
}

TEST(VidermanAlgebraic, SquareOfScaledIdentity) {
  CRSMatrix iI(2, 2);
  iI.row_ptr = {0, 1, 2};
  iI.col_indices = {0, 1};
  iI.values = {Complex(0, 1), Complex(0, 1)};

  CRSMatrix minusI(2, 2);
  minusI.row_ptr = {0, 1, 2};
  minusI.col_indices = {0, 1};
  minusI.values = {Complex(-1, 0), Complex(-1, 0)};

  CompareCRSMatrices(minusI, RunTask(iI, iI));
}

TEST(VidermanAlgebraic, PermutationTimesTranspose) {
  CRSMatrix P(3, 3);
  P.row_ptr = {0, 1, 2, 3};
  P.col_indices = {1, 2, 0};
  P.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix PT(3, 3);
  PT.row_ptr = {0, 1, 2, 3};
  PT.col_indices = {2, 0, 1};
  PT.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix I(3, 3);
  I.row_ptr = {0, 1, 2, 3};
  I.col_indices = {0, 1, 2};
  I.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CompareCRSMatrices(I, RunTask(P, PT));
}

TEST(VidermanStructural, OutputIsValidCRS) {
  CRSMatrix A(4, 3);
  A.row_ptr = {0, 2, 3, 3, 4};
  A.col_indices = {0, 2, 1, 2};
  A.values = {Complex(1, 1), Complex(2, 0), Complex(0, 3), Complex(-1, 1)};

  CRSMatrix B(3, 5);
  B.row_ptr = {0, 2, 3, 5};
  B.col_indices = {0, 3, 2, 1, 4};
  B.values = {Complex(1, 0), Complex(2, 1), Complex(0, 1), Complex(3, 0), Complex(1, -1)};

  CRSMatrix C = RunTask(A, B);
  EXPECT_TRUE(C.IsValid());
  EXPECT_EQ(C.rows, 4);
  EXPECT_EQ(C.cols, 5);
}

TEST(VidermanStructural, ColIndicesSortedInOutput) {
  CRSMatrix A(2, 3);
  A.row_ptr = {0, 3, 3};
  A.col_indices = {0, 1, 2};
  A.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix B(3, 3);
  B.row_ptr = {0, 1, 2, 3};
  B.col_indices = {2, 0, 1};
  B.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix C = RunTask(A, B);
  for (int i = 0; i < C.rows; ++i) {
    for (int j = C.row_ptr[i]; j < C.row_ptr[i + 1] - 1; ++j) {
      EXPECT_LT(C.col_indices[j], C.col_indices[j + 1])
          << "Строка " << i << ": col_indices не отсортированы на позиции " << j;
    }
  }
}

TEST(VidermanStructural, RowPtrStartsAtZero) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(1, 0), Complex(2, 0)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 1};
  B.values = {Complex(3, 0), Complex(4, 0)};

  CRSMatrix C = RunTask(A, B);
  ASSERT_FALSE(C.row_ptr.empty());
  EXPECT_EQ(C.row_ptr[0], 0);
}

TEST(VidermanStructural, RowPtrLastEqualsNNZ) {
  CRSMatrix A(3, 2);
  A.row_ptr = {0, 2, 2, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(1, 1), Complex(2, 0)};

  CRSMatrix B(2, 3);
  B.row_ptr = {0, 2, 3};
  B.col_indices = {0, 2, 1};
  B.values = {Complex(1, 0), Complex(0, 1), Complex(3, 0)};

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(static_cast<size_t>(C.row_ptr[C.rows]), C.values.size());
}

TEST(VidermanStructural, NoExplicitZerosInOutput) {
  CRSMatrix A(2, 2);
  A.row_ptr = {0, 1, 2};
  A.col_indices = {0, 1};
  A.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix B(2, 2);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {1, 0};
  B.values = {Complex(1, 0), Complex(1, 0)};

  CRSMatrix C = RunTask(A, B);
  for (const auto &v : C.values) {
    EXPECT_GT(std::abs(v), 1e-14) << "В values присутствует явный ноль";
  }
}

TEST(VidermanNonTrivial, RectangularWithAccumulation) {
  CRSMatrix A(2, 3);
  A.row_ptr = {0, 2, 3};
  A.col_indices = {0, 2, 1};
  A.values = {Complex(1, 0), Complex(2, 1), Complex(3, 0)};

  CRSMatrix B(3, 4);
  B.row_ptr = {0, 1, 3, 4};
  B.col_indices = {1, 2, 3, 0};
  B.values = {Complex(1, 1), Complex(2, 0), Complex(1, 1), Complex(3, 0)};

  std::vector<std::vector<Complex>> expected = {{Complex(6, 3), Complex(1, 1), Complex(0, 0), Complex(0, 0)},
                                                {Complex(0, 0), Complex(0, 0), Complex(6, 0), Complex(3, 3)}};
  CompareDense(expected, RunTask(A, B));
}

TEST(VidermanNonTrivial, MultipleRowsContributeToSameColumn) {
  CRSMatrix A(3, 2);
  A.row_ptr = {0, 2, 4, 6};
  A.col_indices = {0, 1, 0, 1, 0, 1};
  A.values = {Complex(1, 0), Complex(0, 1), Complex(2, 0), Complex(0, 2), Complex(3, 0), Complex(0, 3)};

  CRSMatrix B(2, 1);
  B.row_ptr = {0, 1, 2};
  B.col_indices = {0, 0};
  B.values = {Complex(1, 1), Complex(1, -1)};

  CRSMatrix expected(3, 1);
  expected.row_ptr = {0, 1, 2, 3};
  expected.col_indices = {0, 0, 0};
  expected.values = {Complex(2, 2), Complex(4, 4), Complex(6, 6)};

  CompareCRSMatrices(expected, RunTask(A, B));
}

TEST(VidermanNonTrivial, CornerElementsOnly5x5) {
  CRSMatrix A(5, 5);
  A.row_ptr = {0, 1, 1, 1, 1, 2};
  A.col_indices = {0, 4};
  A.values = {Complex(1, 0), Complex(0, 1)};

  CRSMatrix B(5, 5);
  B.row_ptr = {0, 1, 1, 1, 1, 2};
  B.col_indices = {0, 4};
  B.values = {Complex(2, 0), Complex(0, -1)};

  CRSMatrix C = RunTask(A, B);
  EXPECT_EQ(C.rows, 5);
  EXPECT_EQ(C.cols, 5);
  EXPECT_EQ(C.NonZeros(), 2u);

  auto D = ToDense(C);
  EXPECT_NEAR(D[0][0].real(), 2.0, 1e-12);
  EXPECT_NEAR(D[0][0].imag(), 0.0, 1e-12);
  EXPECT_NEAR(D[4][4].real(), 1.0, 1e-12);
  EXPECT_NEAR(D[4][4].imag(), 0.0, 1e-12);
}

TEST(VidermanNonTrivial, DenseRowTimesIdentity) {
  CRSMatrix A(1, 4);
  A.row_ptr = {0, 4};
  A.col_indices = {0, 1, 2, 3};
  A.values = {Complex(1, 1), Complex(2, 0), Complex(3, -1), Complex(0, 4)};

  CRSMatrix I(4, 4);
  I.row_ptr = {0, 1, 2, 3, 4};
  I.col_indices = {0, 1, 2, 3};
  I.values = {Complex(1, 0), Complex(1, 0), Complex(1, 0), Complex(1, 0)};

  CRSMatrix C = RunTask(A, I);
  ASSERT_EQ(C.NonZeros(), 4u);
  auto D = ToDense(C);
  EXPECT_NEAR(D[0][0].real(), 1.0, 1e-12);
  EXPECT_NEAR(D[0][0].imag(), 1.0, 1e-12);
  EXPECT_NEAR(D[0][1].real(), 2.0, 1e-12);
  EXPECT_NEAR(D[0][1].imag(), 0.0, 1e-12);
  EXPECT_NEAR(D[0][2].real(), 3.0, 1e-12);
  EXPECT_NEAR(D[0][2].imag(), -1.0, 1e-12);
  EXPECT_NEAR(D[0][3].real(), 0.0, 1e-12);
  EXPECT_NEAR(D[0][3].imag(), 4.0, 1e-12);
}

TEST(VidermanNonTrivial, MatrixSquaredKnownResult) {
  CRSMatrix J(2, 2);
  J.row_ptr = {0, 1, 2};
  J.col_indices = {1, 0};
  J.values = {Complex(1, 0), Complex(-1, 0)};

  CRSMatrix minusI(2, 2);
  minusI.row_ptr = {0, 1, 2};
  minusI.col_indices = {0, 1};
  minusI.values = {Complex(-1, 0), Complex(-1, 0)};

  CompareCRSMatrices(minusI, RunTask(J, J));
}

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
