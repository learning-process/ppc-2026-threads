#include <gtest/gtest.h>

#include <vector>

#include "rychkova_gauss/common/include/common.hpp"
#include "rychkova_gauss/omp/include/ops_omp.hpp"
#include "rychkova_gauss/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace rychkova_gauss {

class RychkovaGaussPerf : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 2500;
  InType input_data_ = Image(kSize_, std::vector<Pixel>(kSize_, Pixel(0, 0, 0)));

  void SetUp() override {
    for (int i = 0; i < kSize_; i++) {
      input_data_[i][0] = Pixel(255, 255, 255);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    bool equal = true;
    for (int j = 0; j < kSize_ && equal; j++) {
      equal &= output_data[j][0] == Pixel(125, 125, 125);
      equal &= output_data[j][1] == Pixel(61, 61, 61);
      for (int i = 2; i < kSize_ && equal; i++) {
        equal &= output_data[j][i] == Pixel(0, 0, 0);
      }
    }
    return equal;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(RychkovaGaussPerf, R) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, RychkovaGaussSEQ, RychkovaGaussOMP>(PPC_SETTINGS_rychkova_gauss);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RychkovaGaussPerf::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, RychkovaGaussPerf, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace rychkova_gauss
