#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <random>
#include <vector>

#include "gusev_d_double_sort_even_odd_batcher/seq/include/ops_seq.hpp"
#include "performance/include/performance.hpp"
#include "util/include/perf_test_util.hpp"

namespace gusev_d_double_sort_even_odd_batcher_task_threads {

class PerfTestBaseSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  bool CheckTestOutputData(OutType &output_data) final {
    auto expected = input_data;
    std::ranges::sort(expected);
    return output_data == expected;
  }

  InType GetTestInputData() final {
    return input_data;
  }

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto start = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [start] {
      const auto now = std::chrono::high_resolution_clock::now();
      const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start).count();
      return static_cast<double>(ns) * 1e-9;
    };
  }

  InType input_data;
};

class RunPerfTestSEQDescending : public PerfTestBaseSEQ {
 protected:
  void SetUp() override {
    input_data.resize(2000);
    for (size_t i = 0; i < input_data.size(); ++i) {
      input_data[i] = static_cast<double>(input_data.size() - i) + (static_cast<double>(i % 7) * 0.01);
    }
  }
};

class RunPerfTestSEQRandom : public PerfTestBaseSEQ {
 protected:
  void SetUp() override {
    input_data.resize(3000);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(-1000000.0, 1000000.0);
    for (double &value : input_data) {
      value = dist(gen);
    }
  }
};

class RunPerfTestSEQNearlySorted : public PerfTestBaseSEQ {
 protected:
  void SetUp() override {
    input_data.resize(3500);
    for (size_t i = 0; i < input_data.size(); ++i) {
      input_data[i] = static_cast<double>(i) * 0.25;
    }

    // Slightly perturb a sorted sequence to emulate realistic near-sorted input.
    for (size_t i = 0; i + 20 < input_data.size(); i += 25) {
      std::swap(input_data[i], input_data[i + 20]);
    }
  }
};

TEST_P(RunPerfTestSEQDescending, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(RunPerfTestSEQRandom, RunPerfModes) {
  ExecuteTest(GetParam());
}

TEST_P(RunPerfTestSEQNearlySorted, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DoubleSortEvenOddBatcherSEQ>(PPC_SETTINGS_gusev_d_double_sort_even_odd_batcher);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestNameDescending = RunPerfTestSEQDescending::CustomPerfTestName;
const auto kPerfTestNameRandom = RunPerfTestSEQRandom::CustomPerfTestName;
const auto kPerfTestNameNearlySorted = RunPerfTestSEQNearlySorted::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTestsDescending, RunPerfTestSEQDescending, kGtestValues, kPerfTestNameDescending);
INSTANTIATE_TEST_SUITE_P(RunModeTestsRandom, RunPerfTestSEQRandom, kGtestValues, kPerfTestNameRandom);
INSTANTIATE_TEST_SUITE_P(RunModeTestsNearlySorted, RunPerfTestSEQNearlySorted, kGtestValues, kPerfTestNameNearlySorted);

}  // namespace

}  // namespace gusev_d_double_sort_even_odd_batcher_task_threads
