#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <libenvpp/detail/environment.hpp>
#include <memory>
#include <random>
#include <string>

#include "gusev_d_double_sort_even_odd_batcher/all/include/common.hpp"
#include "gusev_d_double_sort_even_odd_batcher/all/include/ops_all.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace gusev_d_double_sort_even_odd_batcher_all_task_threads {

namespace {

constexpr size_t kPerfInputSize = 1 << 15;
constexpr int kPerfThreads = 4;

InType GenerateRandomInput(size_t size, uint64_t seed) {
  std::mt19937_64 generator(seed);
  std::uniform_real_distribution<ValueType> distribution(-1.0e6, 1.0e6);

  InType input(size);
  for (ValueType &value : input) {
    value = distribution(generator);
  }

  return input;
}

InType GenerateMixedInput(size_t size) {
  InType input(size);
  for (size_t i = 0; i < size; ++i) {
    const auto sign = (i & 1U) == 0U ? 1.0 : -1.0;
    input[i] = sign * static_cast<ValueType>((size - i) % 2048);
  }

  return input;
}

void SetSteadyPerfTimer(ppc::performance::PerfAttr &perf_attrs) {
  const auto started = std::chrono::high_resolution_clock::now();
  perf_attrs.current_timer = [started] {
    const auto now = std::chrono::high_resolution_clock::now();
    const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - started).count();
    return static_cast<double>(ns) * 1e-9;
  };
}

class AllThreadCountGuard {
 public:
  explicit AllThreadCountGuard(int thread_count) : scoped_("PPC_NUM_THREADS", std::to_string(thread_count)) {}

 private:
  env::detail::set_scoped_environment_variable scoped_;
};

class GusevDoubleSortEvenOddBatcherAllPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    thread_guard_ = std::make_unique<AllThreadCountGuard>(kPerfThreads);
    input_data_ = GenerateRandomInput(kPerfInputSize, 20260320);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto expected = input_data_;
    std::ranges::sort(expected);
    return output_data == expected;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) final {
    SetSteadyPerfTimer(perf_attrs);
  }

 private:
  InType input_data_;
  std::unique_ptr<AllThreadCountGuard> thread_guard_;
};

class GusevDoubleSortEvenOddBatcherAllMixedPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    thread_guard_ = std::make_unique<AllThreadCountGuard>(kPerfThreads);
    input_data_ = GenerateMixedInput(kPerfInputSize);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto expected = input_data_;
    std::ranges::sort(expected);
    return output_data == expected;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) final {
    SetSteadyPerfTimer(perf_attrs);
  }

 private:
  InType input_data_;
  std::unique_ptr<AllThreadCountGuard> thread_guard_;
};

TEST_P(GusevDoubleSortEvenOddBatcherAllPerfTest, RunPerfModesRandom) {
  ExecuteTest(GetParam());
}

TEST_P(GusevDoubleSortEvenOddBatcherAllMixedPerfTest, RunPerfModesMixed) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DoubleSortEvenOddBatcherALL>(PPC_SETTINGS_gusev_d_double_sort_even_odd_batcher);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GusevDoubleSortEvenOddBatcherAllPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GusevDoubleSortEvenOddBatcherAllPerfTest, kGtestValues, kPerfTestName);

INSTANTIATE_TEST_SUITE_P(RunModeTests, GusevDoubleSortEvenOddBatcherAllMixedPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gusev_d_double_sort_even_odd_batcher_all_task_threads
