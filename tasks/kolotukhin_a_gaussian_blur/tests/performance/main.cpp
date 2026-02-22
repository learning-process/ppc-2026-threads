#include <gtest/gtest.h>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "kolotukhin_a_gaussian_blur/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kolotukhin_a_gaussian_blur {

class KolotukhinAGaussinBlurePerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType input_data_;

  void SetUp() override {
    // input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() >= 0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KolotukhinAGaussinBlurePerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KolotukhinAGaussinBlureSEQ>(PPC_SETTINGS_kolotukhin_a_gaussian_blur);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KolotukhinAGaussinBlurePerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KolotukhinAGaussinBlurePerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kolotukhin_a_gaussian_blur
