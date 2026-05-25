#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <utility>

#include "gusev_d_double_sort_even_odd_batcher_omp/common/include/common.hpp"
#include "gusev_d_double_sort_even_odd_batcher_omp/omp/include/ops_omp.hpp"

namespace {

using gusev_d_double_sort_even_odd_batcher_omp_task_threads::DoubleSortEvenOddBatcherOMP;
using gusev_d_double_sort_even_odd_batcher_omp_task_threads::InType;
using gusev_d_double_sort_even_odd_batcher_omp_task_threads::OutType;
using gusev_d_double_sort_even_odd_batcher_omp_task_threads::ValueType;

constexpr size_t kPerfInputSize = 1 << 13;

struct PerfRunResult {
  bool ok = false;
  const char *failed_stage = "";
  OutType output;
  std::chrono::duration<double> elapsed{};

  static PerfRunResult MakeFailure(const char *stage) {
    PerfRunResult result;
    result.failed_stage = stage;
    return result;
  }

  static PerfRunResult MakeSuccess(OutType result_output, std::chrono::duration<double> result_elapsed) {
    PerfRunResult result;
    result.ok = true;
    result.output = std::move(result_output);
    result.elapsed = result_elapsed;
    return result;
  }
};

class GusevDoubleSortEvenOddBatcherOmpEnabledPerf : public ::testing::TestWithParam<int> {};

InType GenerateRandomInput(size_t size, uint64_t seed) {
  std::mt19937_64 generator(seed);
  std::uniform_real_distribution<ValueType> distribution(-1.0e6, 1.0e6);

  InType input(size);
  for (ValueType &value : input) {
    value = distribution(generator);
  }

  return input;
}

InType GenerateDescendingInput(size_t size) {
  InType input(size);
  for (size_t i = 0; i < size; ++i) {
    input[i] = static_cast<ValueType>(size - i);
  }

  return input;
}

InType GenerateNearlySortedInput(size_t size) {
  InType input(size);
  for (size_t i = 0; i < size; ++i) {
    input[i] = static_cast<ValueType>(i);
  }

  for (size_t i = 1; i < size; i += 64) {
    std::swap(input[i - 1], input[i]);
  }

  return input;
}

InType GenerateDuplicateHeavyInput(size_t size) {
  InType input(size);
  for (size_t i = 0; i < size; ++i) {
    input[i] = static_cast<ValueType>((i % 17) - 8);
  }

  return input;
}

PerfRunResult ExecutePerfCase(const InType &input) {
  DoubleSortEvenOddBatcherOMP task(input);
  if (!task.Validation()) {
    return PerfRunResult::MakeFailure("Validation");
  }
  if (!task.PreProcessing()) {
    return PerfRunResult::MakeFailure("PreProcessing");
  }

  const auto started = std::chrono::steady_clock::now();
  if (!task.Run()) {
    return PerfRunResult::MakeFailure("Run");
  }
  const auto finished = std::chrono::steady_clock::now();

  if (!task.PostProcessing()) {
    return PerfRunResult::MakeFailure("PostProcessing");
  }

  return PerfRunResult::MakeSuccess(task.GetOutput(), std::chrono::duration<double>(finished - started));
}

void RunPerfCase(const InType &input) {
  auto expected = input;
  std::ranges::sort(expected);

  const auto result = ExecutePerfCase(input);
  ASSERT_TRUE(result.ok) << result.failed_stage;
  EXPECT_EQ(result.output, expected);
  std::cout << "omp_run_time_sec:" << result.elapsed.count() << '\n';
}

TEST_P(GusevDoubleSortEvenOddBatcherOmpEnabledPerf, RunPerfTestOMPDescending) {
  RunPerfCase(GenerateDescendingInput(kPerfInputSize));
}

TEST_P(GusevDoubleSortEvenOddBatcherOmpEnabledPerf, RunPerfTestOMPRandom) {
  RunPerfCase(GenerateRandomInput(kPerfInputSize, 20260320));
}

TEST_P(GusevDoubleSortEvenOddBatcherOmpEnabledPerf, RunPerfTestOMPNearlySorted) {
  RunPerfCase(GenerateNearlySortedInput(kPerfInputSize));
}

TEST_P(GusevDoubleSortEvenOddBatcherOmpEnabledPerf, RunPerfTestOMPDuplicateHeavy) {
  RunPerfCase(GenerateDuplicateHeavyInput(kPerfInputSize));
}

std::string PrintOmpPerformanceParamName(const ::testing::TestParamInfo<int> &info) {
  static_cast<void>(info);
  return "enabled";
}

INSTANTIATE_TEST_SUITE_P(gusev_d_double_sort_even_odd_batcher_omp_perf, GusevDoubleSortEvenOddBatcherOmpEnabledPerf,
                         ::testing::Values(0), PrintOmpPerformanceParamName);

}  // namespace
