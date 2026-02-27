#include <gtest/gtest.h>

#include <random>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"
#include "kutuzov_i_convex_hull_jarvis/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kutuzov_i_convex_hull_jarvis {

class KutuzovIRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 5000;
  InType input_data_;
  OutType expected_output_;

  void SetUp() override {
    input_data_ = {};

    // NOLINT(cert-msc51-cpp)
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> dist(-1000, 1000);

    for (int i = 0; i < kCount_; i++) {
      double random_x = dist(rng);
      double random_y = dist(rng);

      input_data_.emplace_back(random_x, random_y);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KutuzovIRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KutuzovITestConvexHullSEQ>(PPC_SETTINGS_kutuzov_i_convex_hull_jarvis);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KutuzovIRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KutuzovIRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kutuzov_i_convex_hull_jarvis
