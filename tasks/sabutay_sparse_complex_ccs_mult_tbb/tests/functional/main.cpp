#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "sabutay_sparse_complex_ccs_mult_tbb/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult_tbb/tbb/include/ops_tbb.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

OutType RunTask(const InType &input) {
  SabutayASparseComplexCcsMultTBB task(input);

  EXPECT_EQ(task.GetDynamicTypeOfTask(), ppc::task::TypeOfTask::kTBB);
  EXPECT_TRUE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());
  EXPECT_TRUE(task.Run());
  EXPECT_TRUE(task.PostProcessing());
  EXPECT_TRUE(IsValidCcs(task.GetOutput()));

  return task.GetOutput();
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, MultipliesSquareMatrices) {
  const SparseMatrixCCS lhs = DenseToCcs({
      {Complex{1.0, 2.0}, Complex{}, Complex{}},
      {Complex{}, Complex{-1.0, 1.0}, Complex{}},
      {Complex{}, Complex{}, Complex{2.0, -1.0}},
  });
  const SparseMatrixCCS rhs = DenseToCcs({
      {Complex{2.0, -1.0}, Complex{}, Complex{}},
      {Complex{}, Complex{3.0, 2.0}, Complex{}},
      {Complex{}, Complex{}, Complex{-1.0, 4.0}},
  });

  const SparseMatrixCCS expected = DenseToCcs({
      {Complex{4.0, 3.0}, Complex{}, Complex{}},
      {Complex{}, Complex{-5.0, 1.0}, Complex{}},
      {Complex{}, Complex{}, Complex{2.0, 9.0}},
  });

  const OutType actual = RunTask(std::make_tuple(lhs, rhs));
  EXPECT_TRUE(AreMatricesEqual(actual, expected));
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, MultipliesRectangularMatrices) {
  const SparseMatrixCCS lhs = DenseToCcs({
      {Complex{1.0, 0.0}, Complex{}, Complex{2.0, -1.0}},
      {Complex{}, Complex{3.0, 1.0}, Complex{}},
  });
  const SparseMatrixCCS rhs = DenseToCcs({
      {Complex{1.0, 2.0}, Complex{}},
      {Complex{}, Complex{2.0, 0.0}},
      {Complex{4.0, -1.0}, Complex{1.0, 1.0}},
  });

  const SparseMatrixCCS expected = MultiplyCcsReference(lhs, rhs);

  const OutType actual = RunTask(std::make_tuple(lhs, rhs));
  EXPECT_TRUE(AreMatricesEqual(actual, expected));
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, RemovesZeroEntriesAfterCancellation) {
  const SparseMatrixCCS lhs = DenseToCcs({
      {Complex{1.0, 0.0}, Complex{1.0, 0.0}},
      {Complex{1.0, 0.0}, Complex{-1.0, 0.0}},
  });
  const SparseMatrixCCS rhs = DenseToCcs({
      {Complex{1.0, 0.0}},
      {Complex{-1.0, 0.0}},
  });

  const SparseMatrixCCS expected = DenseToCcs({
      {Complex{}},
      {Complex{2.0, 0.0}},
  });

  const OutType actual = RunTask(std::make_tuple(lhs, rhs));
  EXPECT_TRUE(AreMatricesEqual(actual, expected));
  EXPECT_EQ(actual.row_ind, std::vector<int>({1}));
  ASSERT_EQ(actual.col_ptr.size(), 2U);
  EXPECT_EQ(actual.col_ptr[1], 1);
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, HandlesUnsortedColumnsAndDuplicateRows) {
  SparseMatrixCCS lhs;
  lhs.rows = 3;
  lhs.cols = 2;
  lhs.col_ptr = {0, 2, 5};
  lhs.row_ind = {2, 0, 1, 1, 0};
  lhs.values = {
      Complex{1.0, 1.0}, Complex{2.0, -1.0}, Complex{1.0, 0.0}, Complex{2.0, 0.0}, Complex{-1.0, 2.0},
  };

  SparseMatrixCCS rhs;
  rhs.rows = 2;
  rhs.cols = 2;
  rhs.col_ptr = {0, 2, 4};
  rhs.row_ind = {1, 0, 1, 0};
  rhs.values = {
      Complex{1.0, 0.0},
      Complex{2.0, 0.0},
      Complex{-1.0, 1.0},
      Complex{3.0, -2.0},
  };

  const SparseMatrixCCS expected = MultiplyCcsReference(lhs, rhs);

  const OutType actual = RunTask(std::make_tuple(lhs, rhs));
  EXPECT_TRUE(AreMatricesEqual(actual, expected));
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, SupportsDegenerateZeroInnerDimension) {
  const SparseMatrixCCS lhs = MakeZeroMatrix(3, 0);
  const SparseMatrixCCS rhs = MakeZeroMatrix(0, 4);
  const SparseMatrixCCS expected = MakeZeroMatrix(3, 4);

  const OutType actual = RunTask(std::make_tuple(lhs, rhs));
  EXPECT_TRUE(AreMatricesEqual(actual, expected));
  EXPECT_EQ(actual.rows, 3);
  EXPECT_EQ(actual.cols, 4);
  EXPECT_TRUE(actual.row_ind.empty());
  EXPECT_TRUE(actual.values.empty());
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, ValidationFailsForDimensionMismatchButPipelineIsSafe) {
  const SparseMatrixCCS lhs = BuildDeterministicMatrix(3, 2, 2, 1);
  const SparseMatrixCCS rhs = BuildDeterministicMatrix(3, 2, 2, 5);

  SabutayASparseComplexCcsMultTBB task(std::make_tuple(lhs, rhs));

  EXPECT_FALSE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());
  EXPECT_TRUE(task.Run());
  EXPECT_TRUE(task.PostProcessing());
  EXPECT_TRUE(IsValidCcs(task.GetOutput()));
  EXPECT_EQ(task.GetOutput().rows, 0);
  EXPECT_EQ(task.GetOutput().cols, 0);
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, ValidationFailsForBrokenCcsButPipelineIsSafe) {
  SparseMatrixCCS broken;
  broken.rows = 2;
  broken.cols = 2;
  broken.col_ptr = {0, 2, 1};
  broken.row_ind = {0, 1};
  broken.values = {Complex{1.0, 0.0}, Complex{2.0, 0.0}};

  const SparseMatrixCCS rhs = DenseToCcs({
      {Complex{1.0, 0.0}, Complex{}},
      {Complex{}, Complex{1.0, 0.0}},
  });

  SabutayASparseComplexCcsMultTBB task(std::make_tuple(broken, rhs));

  EXPECT_FALSE(task.Validation());
  EXPECT_TRUE(task.PreProcessing());
  EXPECT_TRUE(task.Run());
  EXPECT_TRUE(task.PostProcessing());
  EXPECT_TRUE(IsValidCcs(task.GetOutput()));
  EXPECT_TRUE(task.GetOutput().row_ind.empty());
  EXPECT_TRUE(task.GetOutput().values.empty());
}

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
