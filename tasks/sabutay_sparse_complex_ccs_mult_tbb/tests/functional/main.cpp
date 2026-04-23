#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "sabutay_sparse_complex_ccs_mult_tbb/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult_tbb/tbb/include/ops_tbb.hpp"
#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

struct TaskExecutionResult {
  ppc::task::TypeOfTask type = ppc::task::TypeOfTask::kTBB;
  bool validation = false;
  bool preprocessing = false;
  bool run = false;
  bool postprocessing = false;
  OutType output;
};

[[nodiscard]] TaskExecutionResult ExecuteTask(const InType &input) {
  SabutayASparseComplexCcsMultTBB task(input);

  TaskExecutionResult result;
  result.type = task.GetDynamicTypeOfTask();
  result.validation = task.Validation();
  result.preprocessing = task.PreProcessing();
  result.run = task.Run();
  result.postprocessing = task.PostProcessing();
  result.output = task.GetOutput();
  return result;
}

[[nodiscard]] bool HasSuccessfulExecution(const TaskExecutionResult &result) {
  return result.type == ppc::task::TypeOfTask::kTBB && result.validation && result.preprocessing && result.run &&
         result.postprocessing && IsValidCcs(result.output);
}

[[nodiscard]] bool HasSafeRejectedExecution(const TaskExecutionResult &result) {
  return result.type == ppc::task::TypeOfTask::kTBB && !result.validation && result.preprocessing && result.run &&
         result.postprocessing && IsValidCcs(result.output);
}

[[nodiscard]] bool IsZeroMatrixOutput(const SparseMatrixCCS &matrix, int rows, int cols) {
  return AreMatricesEqual(matrix, MakeZeroMatrix(rows, cols));
}

[[nodiscard]] bool IsSingleValueColumn(const SparseMatrixCCS &matrix, int row, const Complex &value) {
  return matrix.cols == 1 && matrix.col_ptr == std::vector<int>({0, 1}) && matrix.row_ind == std::vector<int>({row}) &&
         matrix.values.size() == 1U && IsNearZero(matrix.values[0] - value);
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

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSuccessfulExecution(result));
  EXPECT_TRUE(AreMatricesEqual(result.output, expected));
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

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSuccessfulExecution(result));
  EXPECT_TRUE(AreMatricesEqual(result.output, expected));
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

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSuccessfulExecution(result));
  EXPECT_TRUE(AreMatricesEqual(result.output, expected));
  EXPECT_TRUE(IsSingleValueColumn(result.output, 1, Complex{2.0, 0.0}));
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

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSuccessfulExecution(result));
  EXPECT_TRUE(AreMatricesEqual(result.output, expected));
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, SupportsDegenerateZeroInnerDimension) {
  const SparseMatrixCCS lhs = MakeZeroMatrix(3, 0);
  const SparseMatrixCCS rhs = MakeZeroMatrix(0, 4);

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSuccessfulExecution(result));
  EXPECT_TRUE(IsZeroMatrixOutput(result.output, 3, 4));
}

TEST(SabutayASparseComplexCcsMultTBBFunctional, ValidationFailsForDimensionMismatchButPipelineIsSafe) {
  const SparseMatrixCCS lhs = BuildDeterministicMatrix(3, 2, 2, 1);
  const SparseMatrixCCS rhs = BuildDeterministicMatrix(3, 2, 2, 5);

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSafeRejectedExecution(result));
  EXPECT_TRUE(IsZeroMatrixOutput(result.output, 0, 0));
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

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(broken, rhs));
  ASSERT_TRUE(HasSafeRejectedExecution(result));
  EXPECT_TRUE(IsZeroMatrixOutput(result.output, 0, 0));
}

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
