#include <gtest/gtest.h>

#include <cstddef>

#include "morozova_s_strassen_multiplication/common/include/common.hpp"
#include "morozova_s_strassen_multiplication/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace morozova_s_strassen_multiplication {

class MorozovaSStrassenMultiplicationPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int size = 128;
    input_data_ = {static_cast<double>(size)};

    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        input_data_.push_back(static_cast<double>((i * size) + j + 1));
      }
    }

    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        input_data_.push_back(static_cast<double>(((i + j) * 2) + 1));
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int n = static_cast<int>(input_data_[0]);

    if (output_data.empty() || static_cast<int>(output_data[0]) != n) {
      return false;
    }

    double sum = 0.0;
    for (size_t i = 1; i < output_data.size(); ++i) {
      sum += output_data[i];
    }

    return sum > 0.0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(MorozovaSStrassenMultiplicationPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, MorozovaSStrassenMultiplicationSEQ>(
    PPC_SETTINGS_morozova_s_strassen_multiplication);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = MorozovaSStrassenMultiplicationPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(StrassenMultiplicationPerfTests, MorozovaSStrassenMultiplicationPerfTest, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace morozova_s_strassen_multiplication
