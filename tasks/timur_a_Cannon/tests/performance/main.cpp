#include <gtest/gtest.h>

#include <cstddef>

#include "timur_a_Cannon/common/include/common.hpp"
#include "timur_a_Cannon/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace timur_a_Cannon {

class TimurACannonPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const std::size_t kSize_ = 512;
  InType input_data_{};

  void SetUp() override {
    input_data_.A = Matrix(kSize_, 1.0);
    input_data_.B = Matrix(kSize_, 1.0);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.n != kSize_) {
      return false;
    }
    if (output_data.data.size() != kSize_ * kSize_) {
      return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TimurACannonPerfTest, RunPerfSeq) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, TimurACannonSEQ>(PPC_SETTINGS_timur_a_Cannon);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TimurACannonPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TimurACannonPerfTest, kGtestValues, kPerfTestName);

}  // namespace timur_a_Cannon
