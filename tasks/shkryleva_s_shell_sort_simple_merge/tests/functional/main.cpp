#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/seq/include/ops_seq.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

namespace {

TEST(ShkrylevaSShellMergeSEQFunctional, BasicEven) {
  std::vector<int> input{7, 0, 2, 4, 2, 6, 3, 6};
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, NegativeOdd) {
  std::vector<int> input{29, -76, -4, -56, 0, 44, 3};
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, Duplicates) {
  std::vector<int> input{1, 1, 1, 1, 1, 1};
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, Single) {
  std::vector<int> input{222};
  std::vector<int> expected = input;

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, Empty) {
  std::vector<int> input{};
  std::vector<int> expected{};

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, AlreadySorted) {
  std::vector<int> input{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<int> expected = input;

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, ReverseSorted) {
  std::vector<int> input{8, 7, 6, 5, 4, 3, 2, 1};
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

TEST(ShkrylevaSShellMergeSEQFunctional, Small) {
  std::vector<int> input{2, 1};
  std::vector<int> expected = input;
  std::sort(expected.begin(), expected.end());

  auto task = std::make_shared<ShkrylevaSShellMergeSEQ>(input);

  ASSERT_TRUE(task->Validation());
  ASSERT_TRUE(task->PreProcessing());
  ASSERT_TRUE(task->Run());
  ASSERT_TRUE(task->PostProcessing());
  ASSERT_EQ(task->GetOutput(), expected);
}

}  // namespace

}  // namespace shkryleva_s_shell_sort_simple_merge
