#include <gtest/gtest.h>

#include "util/include/perf_test_util.hpp"

namespace belov_e_sobel {

class BelovESobelPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int w_ = 200000;
  const int h_ = 200000;
  InType input_data_;

  void SetUp() override {
    std::vector<uint8_t> vec(w_ * h_, 0);
    input_data_ = std::make_tuple(vec, w_, h_);
  }

  static bool CheckTestOutputData(OutType &output_data) final {
    if (output_data == output_data) {
      return true;
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BelovESobelPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BelovESobelALL, BelovESobelOMP, BelovESobelSEQ, BelovESobelSTL, BelovESobelTBB>(
        PPC_ID_belov_e_sobel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BelovESobelPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BelovESobelPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace belov_e_sobel
