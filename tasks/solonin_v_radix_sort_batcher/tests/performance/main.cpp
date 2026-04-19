#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "solonin_v_radix_sort_batcher/omp/include/ops_omp.hpp"
#include "solonin_v_radix_sort_batcher/seq/include/ops_seq.hpp"
#include "solonin_v_radix_sort_batcher/stl/include/ops_stl.hpp"
#include "solonin_v_radix_sort_batcher/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace solonin_v_radix_sort_batcher {

class RadixBatcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int kSize = 8388608;
    const int kA = 1664525;
    const int kC = 1013904223;
    const int kM = INT32_MAX;
    int seed = 31337;
    input_data_.resize(kSize);
    for (int i = 0; i < kSize; ++i) {
      seed = (kA * seed + kC) % kM;
      input_data_[i] = (seed % 4000) - 2000;
    }
  }

  bool CheckTestOutputData(OutType &output) final {
    for (size_t i = 1; i < output.size(); ++i) {
      if (output[i - 1] > output[i]) return false;
    }
    return true;
  }

  InType GetTestInputData() final { return input_data_; }
};

TEST_P(RadixBatcherPerfTests, RunPerfModes) { ExecuteTest(GetParam()); }

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RadixSortBatcherSEQ, RadixSortBatcherOMP, RadixSortBatcherTBB,
                                RadixSortBatcherSTL>(PPC_SETTINGS_solonin_v_radix_sort_batcher);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = RadixBatcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RadixBatcherPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace solonin_v_radix_sort_batcher
