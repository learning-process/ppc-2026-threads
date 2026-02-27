#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "shekhirev_v_hoare_batcher_sort_seq/common/include/common.hpp"
#include "shekhirev_v_hoare_batcher_sort_seq/seq/include/ops_seq.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

namespace {
std::vector<int> GenerateRandomVector(size_t size) {
  std::vector<int> res(size);
  // NOLINTNEXTLINE(cert-msc51-cpp)
  std::mt19937 gen(73);
  std::uniform_int_distribution<int> dist(-1000, 1000);
  for (size_t i = 0; i < size; ++i) {
    res[i] = dist(gen);
  }
  return res;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void RunFuncTest(const std::vector<int> &in) {
  ShekhirevHoareBatcherSortSEQ task(in);
  ASSERT_TRUE(task.Validation());
  ASSERT_TRUE(task.PreProcessing());
  ASSERT_TRUE(task.Run());
  ASSERT_TRUE(task.PostProcessing());

  std::vector<int> out = task.GetOutput();
  std::vector<int> expected = in;
  std::ranges::sort(expected);

  ASSERT_EQ(out, expected);
}
}  // namespace

TEST(ShekhirevVHoareBatcherSortSEQ, TestEmpty) {
  RunFuncTest({});
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestSingleElement) {
  RunFuncTest({42});
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestSorted) {
  RunFuncTest({1, 2, 3, 4, 5, 6, 7, 8});
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestReverseSorted) {
  RunFuncTest({8, 7, 6, 5, 4, 3, 2, 1});
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestRandomPowerOf2) {
  RunFuncTest(GenerateRandomVector(64));
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestRandomNotPowerOf2) {
  RunFuncTest(GenerateRandomVector(53));
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestLarge) {
  RunFuncTest(GenerateRandomVector(1000));
}

TEST(ShekhirevVHoareBatcherSortSEQ, TestIdenticalElements) {
  RunFuncTest({5, 5, 5, 5, 5, 5, 5, 5});
}

}  // namespace shekhirev_v_hoare_batcher_sort_seq
