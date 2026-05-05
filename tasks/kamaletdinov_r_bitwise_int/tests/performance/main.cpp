#include <gtest/gtest.h>

#include "kamaletdinov_r_bitwise_int/all/include/ops_all.hpp"
#include "kamaletdinov_r_bitwise_int/common/include/common.hpp"
#include "kamaletdinov_r_bitwise_int/omp/include/ops_omp.hpp"
#include "kamaletdinov_r_bitwise_int/seq/include/ops_seq.hpp"
#include "kamaletdinov_r_bitwise_int/stl/include/ops_stl.hpp"
#include "kamaletdinov_r_bitwise_int/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace kamaletdinov_r_bitwise_int {

class KamaletdinovRBitwiseIntRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 1000000;
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

TEST_P(KamaletdinovRBitwiseIntRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KamaletdinovRBitwiseIntALL, KamaletdinovRBitwiseIntOMP,
                                KamaletdinovRBitwiseIntSEQ, KamaletdinovRBitwiseIntSTL, KamaletdinovRBitwiseIntTBB>(
        PPC_SETTINGS_kamaletdinov_r_bitwise_int);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KamaletdinovRBitwiseIntRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KamaletdinovRBitwiseIntRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kamaletdinov_r_bitwise_int
