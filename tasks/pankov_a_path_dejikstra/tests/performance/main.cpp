#include <gtest/gtest.h>

#include "pankov_a_path_dejikstra/all/include/ops_all.hpp"
#include "pankov_a_path_dejikstra/common/include/common.hpp"
#include "pankov_a_path_dejikstra/omp/include/ops_omp.hpp"
#include "pankov_a_path_dejikstra/seq/include/ops_seq.hpp"
#include "pankov_a_path_dejikstra/stl/include/ops_stl.hpp"
#include "pankov_a_path_dejikstra/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace pankov_a_path_dejikstra {

class PankovAPathDejikstraRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
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

TEST_P(PankovAPathDejikstraRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PankovAPathDejikstraALL, PankovAPathDejikstraOMP, PankovAPathDejikstraSEQ,
                                PankovAPathDejikstraSTL, PankovAPathDejikstraTBB>(PPC_SETTINGS_pankov_a_path_dejikstra);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PankovAPathDejikstraRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PankovAPathDejikstraRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace pankov_a_path_dejikstra
