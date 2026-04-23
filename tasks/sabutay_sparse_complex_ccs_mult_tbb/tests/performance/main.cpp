#include <gtest/gtest.h>

#include <tuple>

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

TEST(SabutayASparseComplexCcsMultTBBPerformance, MatchesReferenceOnDeterministicLargeCase) {
  const SparseMatrixCCS lhs = BuildDeterministicMatrix(160, 120, 6, 2);
  const SparseMatrixCCS rhs = BuildDeterministicMatrix(120, 140, 5, 9);
  const SparseMatrixCCS expected = MultiplyCcsReference(lhs, rhs);

  const TaskExecutionResult result = ExecuteTask(std::make_tuple(lhs, rhs));
  ASSERT_TRUE(HasSuccessfulExecution(result));
  EXPECT_TRUE(AreMatricesEqual(result.output, expected));
}

TEST(SabutayASparseComplexCcsMultTBBPerformance, ProducesStableResultAcrossRepeatedRuns) {
  const SparseMatrixCCS lhs = BuildDeterministicMatrix(96, 96, 4, 4);
  const SparseMatrixCCS rhs = BuildDeterministicMatrix(96, 96, 4, 11);

  const TaskExecutionResult first_result = ExecuteTask(std::make_tuple(lhs, rhs));
  const TaskExecutionResult second_result = ExecuteTask(std::make_tuple(lhs, rhs));

  ASSERT_TRUE(HasSuccessfulExecution(first_result));
  ASSERT_TRUE(HasSuccessfulExecution(second_result));
  EXPECT_TRUE(AreMatricesEqual(first_result.output, second_result.output));
}

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
