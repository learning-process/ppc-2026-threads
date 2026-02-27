#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "common/include/common.hpp"
#include "seq/include/ops_seq.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

namespace {
std::vector<int> GenerateRandomVector(size_t size) {
  std::vector<int> res(size);
  std::mt19937 gen(42);
  std::uniform_int_distribution<int> dist(-1000, 1000);
  for (size_t i = 0; i < size; ++i) {
    res[i] = dist(gen);
  }
  return res;
}

void RunFuncTest(const std::vector<int> &in) {
  ShekhirevHoareBatcherSortSEQ task(in);
  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());
  ASSERT_TRUE(task.Run());
  ASSERT_TRUE(task.PostProcessing());

  std::vector<int> out = task.GetOutput();
  std::vector<int> expected = in;
  std::sort(expected.begin(), expected.end());

  ASSERT_EQ(out, expected);
}
}  // namespace

TEST(ShekhirevVHoareBatcherSortSEQ, Test_Empty) {
  RunFuncTest({});
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_SingleElement) {
  RunFuncTest({42});
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_Sorted) {
  RunFuncTest({1, 2, 3, 4, 5, 6, 7, 8});
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_ReverseSorted) {
  RunFuncTest({8, 7, 6, 5, 4, 3, 2, 1});
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_Random_PowerOf2) {
  RunFuncTest(GenerateRandomVector(64));
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_Random_NotPowerOf2) {
  RunFuncTest(GenerateRandomVector(53));
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_Large) {
  RunFuncTest(GenerateRandomVector(1000));
}

TEST(ShekhirevVHoareBatcherSortSEQ, Test_IdenticalElements) {
  RunFuncTest({5, 5, 5, 5, 5, 5, 5, 5});
}
}  // namespace shekhirev_v_hoare_batcher_sort_seq
