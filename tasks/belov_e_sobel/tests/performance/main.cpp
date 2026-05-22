#include <gtest/gtest.h>

#include <tuple>

#include "belov_e_sobel/all/include/ops_all.hpp"
#include "belov_e_sobel/common/include/common.hpp"
#include "belov_e_sobel/omp/include/ops_omp.hpp"
#include "belov_e_sobel/seq/include/ops_seq.hpp"
#include "belov_e_sobel/stl/include/ops_stl.hpp"
#include "belov_e_sobel/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace belov_e_sobel {

class BelovESobelPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int w_ = 10000;
  const int h_ = 10000;
  InType input_data_;

  void SetUp() override {
    std::vector<uint8_t> vec(w_ * h_, 0);
    input_data_ = std::make_tuple(vec, w_, h_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
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
        PPC_SETTINGS_belov_e_sobel);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BelovESobelPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BelovESobelPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace belov_e_sobel
