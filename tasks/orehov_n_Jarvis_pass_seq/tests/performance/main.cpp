#include <gtest/gtest.h>
#include <cmath>

#include "orehov_n_Jarvis_pass_seq/common/include/common.hpp"
#include "orehov_n_Jarvis_pass_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace orehov_n_Jarvis_pass_seq {

class OrehovNJarvisPassSEQPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 2000;
  InType input_data_{};

  void SetUp() override {
    double radius = 2000;
    for (int i = 0; i < kCount_ - 1; i++){
      double angle = 2 * 3.14 * i / kCount_;
        
      double x = radius * std::cos(angle);
      double y = radius * std::sin(angle);
        
      input_data_.emplace_back(x, y);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {

    return (output_data.size() > 0);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(OrehovNJarvisPassSEQPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, OrehovNJarvisPassSEQ>(PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = OrehovNJarvisPassSEQPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, OrehovNJarvisPassSEQPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace orehov_n_Jarvis_pass_seq
