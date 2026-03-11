#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

#include "gusev_d_double_sort_even_odd_batcher/seq/include/ops_seq.hpp"

namespace gusev_d_double_sort_even_odd_batcher_task_threads {
namespace {

bool RunTaskPipeline(DoubleSortEvenOddBatcherSEQ &task) {
  return task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
}

OutType RunTask(const InType &input) {
  DoubleSortEvenOddBatcherSEQ task(input);
  if (!RunTaskPipeline(task)) {
    throw std::runtime_error("Task pipeline failed");
  }
  return task.GetOutput();
}

void CheckSortedPermutation(const InType &input, const OutType &output) {
  ASSERT_EQ(output.size(), input.size());
  EXPECT_TRUE(std::ranges::is_sorted(output));
  EXPECT_TRUE(std::ranges::is_permutation(output, input));
}

void CheckMatchesStdSort(const InType &input, const OutType &output) {
  auto expected = input;
  std::ranges::sort(expected);
  EXPECT_EQ(output, expected);
}

InType GenerateRandomInput(size_t size, uint64_t seed) {
  std::mt19937_64 gen(seed);
  std::uniform_real_distribution<double> dist(-1000000.0, 1000000.0);
  InType data(size);
  for (double &value : data) {
    value = dist(gen);
  }
  return data;
}

template <typename Callback>
bool ThrowsRuntimeError(Callback callback) {
  try {
    callback();
  } catch (const std::runtime_error &) {
    return true;
  }
  return false;
}

bool ValidationFailsWhenOutputPrepared() {
  const InType input{3.0, 1.0, 2.0};
  DoubleSortEvenOddBatcherSEQ task(input);
  task.GetOutput() = {42.0};
  return !task.Validation() && ThrowsRuntimeError([&task] { task.Run(); });
}

bool PreProcessingBeforeValidationThrows() {
  const InType input{1.0, 0.0};
  DoubleSortEvenOddBatcherSEQ task(input);
  return ThrowsRuntimeError([&task] { task.PreProcessing(); });
}

bool RunBeforePreProcessingThrows() {
  const InType input{1.0, 0.0};
  DoubleSortEvenOddBatcherSEQ task(input);
  if (!task.Validation()) {
    return false;
  }
  return ThrowsRuntimeError([&task] { task.Run(); });
}

bool PostProcessingBeforeRunThrows() {
  const InType input{1.0, 0.0};
  DoubleSortEvenOddBatcherSEQ task(input);
  if (!task.Validation() || !task.PreProcessing()) {
    return false;
  }
  return ThrowsRuntimeError([&task] { task.PostProcessing(); });
}

bool RunTaskTwiceBeforePostProcessing(const InType &input, OutType &output) {
  DoubleSortEvenOddBatcherSEQ task(input);
  if (!task.Validation() || !task.PreProcessing() || !task.Run() || !task.Run() || !task.PostProcessing()) {
    return false;
  }
  output = task.GetOutput();
  return true;
}

bool RunTaskWithInputSnapshot(const InType &original_input, const InType &replacement_input, OutType &output) {
  DoubleSortEvenOddBatcherSEQ task(original_input);
  if (!task.Validation() || !task.PreProcessing()) {
    return false;
  }

  task.GetInput() = replacement_input;

  if (!task.Run() || !task.PostProcessing()) {
    return false;
  }

  output = task.GetOutput();
  return true;
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, HandlesEmptyVector) {
  const InType input{};
  const auto output = RunTask(input);
  EXPECT_TRUE(output.empty());
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsSingleElement) {
  const InType input{42.5};
  const auto output = RunTask(input);
  EXPECT_EQ(output, input);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsMixedValues) {
  const InType input{3.14, -10.0, 2.71, 0.0, -0.25, 100.1, -99.9};
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsWithDuplicates) {
  const InType input{5.0, 2.0, 5.0, -1.0, 2.0, -1.0, 0.0};
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsReversedSequence) {
  InType input(50);
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = static_cast<double>(input.size() - i);
  }
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsAlreadySortedSequence) {
  InType input(128);
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = static_cast<double>(i) * 0.5;
  }
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsTwoElements) {
  const InType input{10.0, -10.0};
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsAllEqualValues) {
  const InType input(200, 7.25);
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsOddSizeDeterministicRandom) {
  const InType input = GenerateRandomInput(257, 42);
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsEvenSizeDeterministicRandom) {
  const InType input = GenerateRandomInput(512, 1337);
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsLargeDeterministicRandom) {
  const InType input = GenerateRandomInput(4096, 20260308);
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsExtremeValues) {
  const InType input{
      std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(),
      std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::min(),
      -std::numeric_limits<double>::min(),
      0.0,
      -0.0,
      1.0,
      -1.0,
  };
  const auto output = RunTask(input);
  CheckSortedPermutation(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, ValidationFailsWhenOutputAlreadyPrepared) {
  EXPECT_TRUE(ValidationFailsWhenOutputPrepared());
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, ThrowsIfPreProcessingCalledBeforeValidation) {
  EXPECT_TRUE(PreProcessingBeforeValidationThrows());
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, ThrowsIfRunCalledBeforePreProcessing) {
  EXPECT_TRUE(RunBeforePreProcessingThrows());
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, ThrowsIfPostProcessingCalledBeforeRun) {
  EXPECT_TRUE(PostProcessingBeforeRunThrows());
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, AllowsRepeatedRunBeforePostProcessing) {
  const InType input = GenerateRandomInput(333, 20260310);
  const auto reference_output = RunTask(input);

  OutType repeated_run_output;
  ASSERT_TRUE(RunTaskTwiceBeforePostProcessing(input, repeated_run_output));
  CheckSortedPermutation(input, repeated_run_output);
  EXPECT_EQ(repeated_run_output, reference_output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, UsesInputSnapshotFromPreProcessing) {
  const InType original_input{9.0, -1.0, 4.0, 2.0, -5.0, 3.0};
  OutType output;
  ASSERT_TRUE(RunTaskWithInputSnapshot(original_input, {-1000.0, -2000.0, -3000.0}, output));
  CheckMatchesStdSort(original_input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsPrimeSizeDeterministicRandomAndMatchesStdSort) {
  const InType input = GenerateRandomInput(997, 8675309);
  const auto output = RunTask(input);
  CheckMatchesStdSort(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsInterleavedLargeAndTinyMagnitudes) {
  const InType input{
      1.0e308,          -1.0e308,         1.0e-308, -1.0e-308, 5.0e307, -5.0e307, 7.5e-309, -7.5e-309,
      std::numbers::pi, -std::numbers::e, 0.0,      -0.0,      42.0,    -42.0,    9.0e307,  -9.0e307,
  };
  const auto output = RunTask(input);
  CheckMatchesStdSort(input, output);
}

TEST(GusevDoubleSortEvenOddBatcherSEQ, SortsDenseDuplicatePatternAndMatchesStdSort) {
  InType input;
  input.reserve(1200);
  for (size_t i = 0; i < 300; ++i) {
    input.push_back(1.0);
    input.push_back(-1.0);
    input.push_back(0.0);
    input.push_back(1.0);
  }

  const auto output = RunTask(input);
  CheckMatchesStdSort(input, output);
}

}  // namespace
}  // namespace gusev_d_double_sort_even_odd_batcher_task_threads
