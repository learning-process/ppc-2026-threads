#include <gtest/gtest.h>

#include "fatehov_k_gaussian/common/include/common.hpp"
#include "fatehov_k_gaussian/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace fatehov_k_gaussian {

class FatehovKGaussianPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 256;
  InType input_data_;

  void SetUp() override {
    input_data_.image = Image(kCount_, kCount_, 3);
    std::fill(input_data_.image.data.begin(), input_data_.image.data.end(), 128);
    input_data_.sigma = 1.0f;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(FatehovKGaussianPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, FatehovKGaussianSEQ>(PPC_SETTINGS_fatehov_k_gaussian);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(FatehovKPerfTests, FatehovKGaussianPerfTests, kGtestValues,
                         FatehovKGaussianPerfTests::CustomPerfTestName);
}  // namespace
}  // namespace fatehov_k_gaussian
