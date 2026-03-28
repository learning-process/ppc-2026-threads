#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <random>
#include <utility>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"
#include "smyshlaev_a_sle_cg_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace smyshlaev_a_sle_cg_seq {

class SmyshlaevASleCgPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kSystemSize = 512;
  InType input_data_{};
  OutType expected_x_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    CGMatrix a(kSystemSize, CGVector(kSystemSize, 0.0));
    expected_x_.assign(kSystemSize, 1.0);
    CGVector b(kSystemSize, 0.0);

    for (int i = 0; i < kSystemSize; ++i) {
      for (int j = i + 1; j < kSystemSize; ++j) {
        double val = dis(gen);
        a[i][j] = val;
        a[j][i] = val;
      }
    }

    for (int i = 0; i < kSystemSize; ++i) {
      double row_sum = 0.0;
      for (int j = 0; j < kSystemSize; ++j) {
        if (i != j) {
          row_sum += std::abs(a[i][j]);
        }
      }
      a[i][i] = row_sum + 1.0;
    }

    for (int i = 0; i < kSystemSize; ++i) {
      for (int j = 0; j < kSystemSize; ++j) {
        b[i] += a[i][j] * expected_x_[j];
      }
    }

    input_data_.A = std::move(a);
    input_data_.b = std::move(b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_x_.size()) {
      return false;
    }
    const double epsilon = 1e-6;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::isnan(output_data[i]) || std::fabs(output_data[i] - expected_x_[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SmyshlaevASleCgPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SmyshlaevASleCgTaskSEQ>(PPC_SETTINGS_smyshlaev_a_sle_cg_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SmyshlaevASleCgPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SmyshlaevASleCgPerfTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace smyshlaev_a_sle_cg_seq
