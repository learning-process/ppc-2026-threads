#include <gtest/gtest.h>

#include <tuple>

#include "karpich_i_bitwise_batcher/all/include/ops_all.hpp"
#include "karpich_i_bitwise_batcher/common/include/common.hpp"
#include "karpich_i_bitwise_batcher/omp/include/ops_omp.hpp"
#include "karpich_i_bitwise_batcher/seq/include/ops_seq.hpp"
#include "karpich_i_bitwise_batcher/stl/include/ops_stl.hpp"
#include "karpich_i_bitwise_batcher/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace karpich_i_bitwise_batcher {

class KarpichIBitwiseBatcherPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1500000;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KarpichIBitwiseBatcherPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = std::tuple_cat(
    ppc::util::MakeAllPerfTasks<InType, KarpichIBitwiseBatcherSEQ>(PPC_SETTINGS_karpich_i_bitwise_batcher),
    ppc::util::MakeAllPerfTasks<InType, KarpichIBitwiseBatcherSTL>(PPC_SETTINGS_karpich_i_bitwise_batcher),
    ppc::util::MakeAllPerfTasks<InType, KarpichIBitwiseBatcherOMP>(PPC_SETTINGS_karpich_i_bitwise_batcher),
    ppc::util::MakeAllPerfTasks<InType, KarpichIBitwiseBatcherTBB>(PPC_SETTINGS_karpich_i_bitwise_batcher),
    ppc::util::MakeAllPerfTasks<InType, KarpichIBitwiseBatcherALL>(PPC_SETTINGS_karpich_i_bitwise_batcher));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KarpichIBitwiseBatcherPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KarpichIBitwiseBatcherPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace karpich_i_bitwise_batcher
