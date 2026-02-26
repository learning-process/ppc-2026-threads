#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <ranges>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

namespace {

void RunTask(const std::vector<int> &input, const std::vector<int> &expected) {
  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);
  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

void RunAndCheck(const std::vector<int> &input) {
  std::vector<int> expected = input;
  std::ranges::sort(expected.begin(), expected.end());
  RunTask(input, expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, BasicEven) {
  RunAndCheck({7, 0, 2, 4, 2, 6, 3, 6});
}

TEST(ShkrylevaSShellMergeSEQFunctional, NegativeOdd) {
  RunAndCheck({29, -76, -4, -56, 0, 44, 3});
}

TEST(ShkrylevaSShellMergeSEQFunctional, Duplicates) {
  RunAndCheck({1, 1, 1, 1, 1, 1});
}

TEST(ShkrylevaSShellMergeSEQFunctional, Single) {
  RunAndCheck({222});
}

TEST(ShkrylevaSShellMergeSEQFunctional, Empty) {
  RunAndCheck({});
}

TEST(ShkrylevaSShellMergeSEQFunctional, AlreadySorted) {
  RunAndCheck({1, 2, 3, 4, 5, 6, 7, 8});
}

TEST(ShkrylevaSShellMergeSEQFunctional, ReverseSorted) {
  RunAndCheck({8, 7, 6, 5, 4, 3, 2, 1});
}

TEST(ShkrylevaSShellMergeSEQFunctional, Small) {
  RunAndCheck({2, 1});
}
}  // namespace

// namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
