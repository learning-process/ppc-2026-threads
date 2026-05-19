#include <gtest/gtest.h>

#include <cstddef>
#include <utility>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"
#include "timofeev_n_radix_batcher_sort/omp/include/ops_omp.hpp"
#include "timofeev_n_radix_batcher_sort/seq/include/ops_seq.hpp"
#include "timofeev_n_radix_batcher_sort/stl/include/ops_stl.hpp"
#include "timofeev_n_radix_batcher_sort/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

class TimofeevRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    int size = 10000;
    input_data_.resize(size, 0);
    for (int i = 0; std::cmp_less(i, static_cast<int>(input_data_.size())); i++) {
      input_data_[i] = i * (i % 2 == 0 ? 1 : -1);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    bool is_true = true;
    for (size_t i = 0; i < output_data.size() - 1; i++) {
      if (output_data[i] > output_data[i + 1]) {
        is_true = false;
      }
    }
    return is_true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TimofeevRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TimofeevNRadixBatcherSEQ, TimofeevNRadixBatcherTBB, TimofeevNRadixBatcherOMP,
                                TimofeevNRadixBatcherSTL>(PPC_SETTINGS_timofeev_n_radix_batcher_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TimofeevRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TimofeevRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace timofeev_n_radix_batcher_sort_threads
