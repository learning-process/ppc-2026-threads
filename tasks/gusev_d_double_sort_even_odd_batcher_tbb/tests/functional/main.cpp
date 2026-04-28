#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <libenvpp/detail/environment.hpp>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

#include "gusev_d_double_sort_even_odd_batcher_tbb/common/include/common.hpp"
#include "gusev_d_double_sort_even_odd_batcher_tbb/tbb/include/ops_tbb.hpp"

namespace {

using gusev_d_double_sort_even_odd_batcher_tbb_task_threads::DoubleSortEvenOddBatcherTBB;
using gusev_d_double_sort_even_odd_batcher_tbb_task_threads::InType;
using gusev_d_double_sort_even_odd_batcher_tbb_task_threads::OutType;
using gusev_d_double_sort_even_odd_batcher_tbb_task_threads::ValueType;

InType GenerateRandomInput(size_t size, uint64_t seed) {
  std::mt19937_64 generator(seed);
  std::uniform_real_distribution<ValueType> distribution(-1.0e6, 1.0e6);

  InType input(size);
  for (ValueType &value : input) {
    value = distribution(generator);
  }

  return input;
}

InType MakeDenseDuplicateInput() {
  return {5.0, 3.0, 5.0, -1.0, -1.0, 0.0, 3.0, 3.0, 8.0, 8.0, 8.0, -4.0, -4.0, 2.0, 2.0};
}

InType MakeExtremeInput() {
  return {std::numbers::pi,
          -std::numbers::e,
          0.0,
          -0.0,
          42.0,
          -42.0,
          std::numeric_limits<ValueType>::max(),
          std::numeric_limits<ValueType>::lowest(),
          std::numeric_limits<ValueType>::min(),
          -std::numeric_limits<ValueType>::min(),
          1.0e-12,
          -1.0e-12,
          7.5,
          7.5,
          -7.5};
}

InType MakeAlternatingMagnitudeInput() {
  return {1.0e-12, -1.0e12, 2.0e-12, -2.0e12,  3.0e-12, -3.0e12,  4.0e-12, -4.0e12,
          5.0e-12, -5.0e12, 6.0e-12, -6.0e12,  7.0e-12, -7.0e12,  8.0e-12, -8.0e12,
          9.0e-12, -9.0e12, 1.0e13,  -1.0e-13, 2.0e13,  -2.0e-13, 3.0e13,  -3.0e-13};
}

class GusevDoubleSortEvenOddBatcherTbbEnabled : public ::testing::TestWithParam<int> {};

class TbbThreadCountGuard {
 public:
  explicit TbbThreadCountGuard(int thread_count) : scoped_("PPC_NUM_THREADS", std::to_string(thread_count)) {}

 private:
  env::detail::set_scoped_environment_variable scoped_;
};

OutType RunTaskPipeline(const InType &input) {
  DoubleSortEvenOddBatcherTBB task(input);
  if (!task.Validation()) {
    throw std::runtime_error("Validation failed");
  }
  if (!task.PreProcessing()) {
    throw std::runtime_error("PreProcessing failed");
  }
  if (!task.Run()) {
    throw std::runtime_error("Run failed");
  }
  if (!task.PostProcessing()) {
    throw std::runtime_error("PostProcessing failed");
  }

  return task.GetOutput();
}

void CheckMatchesStdSort(const InType &input) {
  auto expected = input;
  std::ranges::sort(expected);

  const auto output = RunTaskPipeline(input);
  EXPECT_EQ(output, expected);
}

bool ValidationRejectsPreparedOutputImpl() {
  DoubleSortEvenOddBatcherTBB task({3.0, 2.0, 1.0});
  task.GetOutput() = {0.0};

  if (task.Validation()) {
    return false;
  }

  try {
    static_cast<void>(task.Run());
  } catch (const std::runtime_error &) {
    return true;
  }

  return false;
}

bool ThrowsIfPreProcessingBeforeValidationImpl() {
  DoubleSortEvenOddBatcherTBB task({1.0});

  try {
    static_cast<void>(task.PreProcessing());
  } catch (const std::runtime_error &) {
    return true;
  }

  return false;
}

bool ThrowsIfRunBeforePreProcessingImpl() {
  DoubleSortEvenOddBatcherTBB task({2.0, 1.0});
  if (!task.Validation()) {
    return false;
  }

  try {
    static_cast<void>(task.Run());
  } catch (const std::runtime_error &) {
    return true;
  }

  return false;
}

bool ThrowsIfPostProcessingBeforeRunImpl() {
  DoubleSortEvenOddBatcherTBB task({2.0, 1.0});
  if (!task.Validation() || !task.PreProcessing()) {
    return false;
  }

  try {
    static_cast<void>(task.PostProcessing());
  } catch (const std::runtime_error &) {
    return true;
  }

  return false;
}

OutType RunTaskTwiceBeforePostProcessing(const InType &input) {
  DoubleSortEvenOddBatcherTBB task(input);
  if (!task.Validation() || !task.PreProcessing() || !task.Run() || !task.Run() || !task.PostProcessing()) {
    throw std::runtime_error("Repeated run pipeline failed");
  }

  return task.GetOutput();
}

OutType RunWithInputMutationAfterPreProcessing(const InType &input, ValueType first_value, ValueType second_value) {
  DoubleSortEvenOddBatcherTBB task(input);
  if (!task.Validation() || !task.PreProcessing()) {
    throw std::runtime_error("Preprocessing pipeline failed");
  }

  auto &input_ref = task.GetInput();
  input_ref[0] = first_value;
  input_ref[1] = second_value;

  if (!task.Run() || !task.PostProcessing()) {
    throw std::runtime_error("Run pipeline failed");
  }

  return task.GetOutput();
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsEmptyInput) {
  CheckMatchesStdSort({});
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsSingleElement) {
  CheckMatchesStdSort({42.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsAlreadySortedInput) {
  CheckMatchesStdSort({-7.0, -2.0, -0.0, 0.0, 1.5, 3.0, 4.0, 9.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsReverseSortedInput) {
  CheckMatchesStdSort({9.0, 4.0, 3.0, 1.5, 0.0, -0.0, -2.0, -7.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsOddSizedInput) {
  CheckMatchesStdSort({3.0, -1.0, 2.0, 0.0, 5.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsDenseDuplicateInput) {
  CheckMatchesStdSort(MakeDenseDuplicateInput());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, SortsAllEqualInput) {
  CheckMatchesStdSort(InType(257, 3.5));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForPrimeSizedRandomInput) {
  CheckMatchesStdSort(GenerateRandomInput(997, 20260320));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForLargeRandomInput) {
  CheckMatchesStdSort(GenerateRandomInput(4096, 20260321));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForExtremesAndSignedZeros) {
  CheckMatchesStdSort(MakeExtremeInput());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForAlternatingMagnitudeInput) {
  CheckMatchesStdSort(MakeAlternatingMagnitudeInput());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortWhenInputSizeIsLessThanThreadCount) {
  const TbbThreadCountGuard guard(8);
  CheckMatchesStdSort({7.0, -4.0, 2.5});
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForOddNumberOfBlocks) {
  const TbbThreadCountGuard guard(5);
  CheckMatchesStdSort(GenerateRandomInput(23, 20260322));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForParallelMergeWithOddCarryBlock) {
  const TbbThreadCountGuard guard(65);
  CheckMatchesStdSort(GenerateRandomInput(65, 20260325));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, MatchesStdSortForSmallRandomInput) {
  CheckMatchesStdSort(GenerateRandomInput(23, 20260323));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, ValidationRejectsPreparedOutput) {
  EXPECT_TRUE(ValidationRejectsPreparedOutputImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, ThrowsIfPreProcessingBeforeValidation) {
  EXPECT_TRUE(ThrowsIfPreProcessingBeforeValidationImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, ThrowsIfRunBeforePreProcessing) {
  EXPECT_TRUE(ThrowsIfRunBeforePreProcessingImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, ThrowsIfPostProcessingBeforeRun) {
  EXPECT_TRUE(ThrowsIfPostProcessingBeforeRunImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, AllowsRepeatedRunBeforePostProcessing) {
  const InType input{9.0, -1.0, 5.0, 3.0, -7.0, 11.0, 0.0, 2.0};
  const auto expected_output = RunTaskPipeline(input);
  EXPECT_EQ(RunTaskTwiceBeforePostProcessing(input), expected_output);
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, UsesInputSnapshotFromPreProcessing) {
  EXPECT_EQ(RunWithInputMutationAfterPreProcessing({5.0, 4.0, 3.0, 2.0, 1.0}, -100.0, 200.0),
            (OutType{1.0, 2.0, 3.0, 4.0, 5.0}));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, KeepsElementMultiplicity) {
  const auto input = MakeDenseDuplicateInput();
  const auto output = RunTaskPipeline(input);

  EXPECT_TRUE(std::ranges::is_sorted(output));
  EXPECT_TRUE(std::ranges::is_permutation(output, input));
}

TEST_P(GusevDoubleSortEvenOddBatcherTbbEnabled, KeepsOutputEmptyAfterRunningOnEmptyInput) {
  const auto output = RunTaskPipeline({});
  EXPECT_TRUE(output.empty());
}

std::string PrintTbbFunctionalParamName(const ::testing::TestParamInfo<int> &info) {
  static_cast<void>(info);
  return "enabled";
}

INSTANTIATE_TEST_SUITE_P(gusev_d_double_sort_even_odd_batcher_tbb_enabled, GusevDoubleSortEvenOddBatcherTbbEnabled,
                         ::testing::Values(0), PrintTbbFunctionalParamName);

}  // namespace
