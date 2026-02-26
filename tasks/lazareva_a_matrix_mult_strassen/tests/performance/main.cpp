#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"
#include "lazareva_a_matrix_mult_strassen/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace lazareva_a_matrix_mult_strassen {

class LazarevaARunPerfTestThreads
    : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kN_ = 512;
  InType input_data_{};
  OutType expected_output_;

  void SetUp() override {
    const int size = kN_ * kN_;
    std::vector<double> a(static_cast<size_t>(size));
    std::vector<double> b(static_cast<size_t>(size));

    for (int i = 0; i < size; ++i) {
      a[static_cast<size_t>(i)] =
          static_cast<double>((i % 7) + 1);
      b[static_cast<size_t>(i)] =
          static_cast<double>(((i * 3 + 5) % 11) + 1);
    }

    input_data_ = MatrixInput{.a = a, .b = b, .n = kN_};

    expected_output_.assign(static_cast<size_t>(size), 0.0);
    for (int row = 0; row < kN_; ++row) {
      for (int k = 0; k < kN_; ++k) {
        for (int col = 0; col < kN_; ++col) {
          expected_output_[static_cast<size_t>((row * kN_) + col)] +=
              a[static_cast<size_t>((row * kN_) + k)] *
              b[static_cast<size_t>((k * kN_) + col)];
        }
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    constexpr double kEps = 1e-6;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_output_[i]) > kEps) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final { return input_data_; }
};

TEST_P(LazarevaARunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, LazarevaATestTaskSEQ>(
        PPC_SETTINGS_lazareva_a_matrix_mult_strassen);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LazarevaARunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, LazarevaARunPerfTestThreads,
                         kGtestValues, kPerfTestName);

}  // namespace

}  // namespace lazareva_a_matrix_mult_strassen