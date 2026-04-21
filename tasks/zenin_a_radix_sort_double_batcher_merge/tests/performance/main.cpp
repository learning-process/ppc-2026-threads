#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

// #include "zenin_a_radix_sort_double_batcher_merge/all/include/ops_all.hpp"
#include "util/include/perf_test_util.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/omp/include/ops_omp.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/seq/include/ops_seq.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/stl/include/ops_stl.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/tbb/include/ops_tbb.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

class ZeninARadixSortDoubleBatcherMergePerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1000000;
  InType input_data_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1e6, 1e6);

    input_data_.resize(kCount_);
    for (int i = 0; i < kCount_; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected_data = input_data_;
    std::ranges::sort(expected_data);
    if (output_data.size() != expected_data.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i] != expected_data[i]) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ZeninARadixSortDoubleBatcherMergePerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZeninARadixSortDoubleBatcherMergeOMP, ZeninARadixSortDoubleBatcherMergeSeqseq,
                                ZeninARadixSortDoubleBatcherMergeTBB, ZeninARadixSortDoubleBatcherMergeSTL>(
        PPC_SETTINGS_zenin_a_radix_sort_double_batcher_merge);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZeninARadixSortDoubleBatcherMergePerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZeninARadixSortDoubleBatcherMergePerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zenin_a_radix_sort_double_batcher_merge
