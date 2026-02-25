#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "sakharov_a_shell_sorting_with_merging_butcher/seq/include/ops_seq.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

namespace {

void RunAndCheck(const std::vector<int> &input) {
  auto task = std::make_shared<SakharovAShellButcherSEQ>(input);
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(SakharovAShellButcherSEQFunctional, BasicEven) {
  RunAndCheck({5, 1, 8, 3, 2, 7, 4, 6});
}

TEST(SakharovAShellButcherSEQFunctional, NegativeOdd) {
  RunAndCheck({9, -2, 4, -7, 0, 3, -1});
}

TEST(SakharovAShellButcherSEQFunctional, Duplicates) {
  RunAndCheck({1, 1, 1, 1, 1, 1});
}

TEST(SakharovAShellButcherSEQFunctional, Single) {
  RunAndCheck({42});
}

TEST(SakharovAShellButcherSEQFunctional, Empty) {
  RunAndCheck({});
}

}  // namespace

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
