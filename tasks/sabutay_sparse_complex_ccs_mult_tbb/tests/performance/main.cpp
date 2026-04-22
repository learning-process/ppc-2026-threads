#include <gtest/gtest.h>

#include <tuple>

#include "sabutay_sparse_complex_ccs_mult_tbb/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_mult_tbb/tbb/include/ops_tbb.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

TEST(SabutayASparseComplexCcsMultTBBPerformance, MatchesReferenceOnDeterministicLargeCase) {
  const SparseMatrixCCS lhs = BuildDeterministicMatrix(160, 120, 6, 2);
  const SparseMatrixCCS rhs = BuildDeterministicMatrix(120, 140, 5, 9);
  const SparseMatrixCCS expected = MultiplyCcsReference(lhs, rhs);

  SabutayASparseComplexCcsMultTBB task(std::make_tuple(lhs, rhs));

  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());
  ASSERT_TRUE(task.Run());
  ASSERT_TRUE(task.PostProcessing());
  ASSERT_TRUE(IsValidCcs(task.GetOutput()));
  EXPECT_TRUE(AreMatricesEqual(task.GetOutput(), expected));
}

TEST(SabutayASparseComplexCcsMultTBBPerformance, ProducesStableResultAcrossRepeatedRuns) {
  const SparseMatrixCCS lhs = BuildDeterministicMatrix(96, 96, 4, 4);
  const SparseMatrixCCS rhs = BuildDeterministicMatrix(96, 96, 4, 11);

  SabutayASparseComplexCcsMultTBB first_task(std::make_tuple(lhs, rhs));
  SabutayASparseComplexCcsMultTBB second_task(std::make_tuple(lhs, rhs));

  ASSERT_TRUE(first_task.Validation());
  ASSERT_TRUE(first_task.PreProcessing());
  ASSERT_TRUE(first_task.Run());
  ASSERT_TRUE(first_task.PostProcessing());

  ASSERT_TRUE(second_task.Validation());
  ASSERT_TRUE(second_task.PreProcessing());
  ASSERT_TRUE(second_task.Run());
  ASSERT_TRUE(second_task.PostProcessing());

  EXPECT_TRUE(AreMatricesEqual(first_task.GetOutput(), second_task.GetOutput()));
}

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
