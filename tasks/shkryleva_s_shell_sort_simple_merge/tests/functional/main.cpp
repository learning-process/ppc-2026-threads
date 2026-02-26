#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

namespace {

void RunAndCheck(const std::vector<int> &input) {
  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
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
