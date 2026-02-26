#include <gtest/gtest.h>

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"
#include "tabalaev_a_matrix_mul_strassen/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace tabalaev_a_matrix_mul_strassen {

class TabalaevAMatrixMulStrassenPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    const size_t kRc = 512;
    const size_t size = kRc * kRc;

    input_data_.a_rows = static_cast<int>(kRc);
    input_data_.a_cols_b_rows = static_cast<int>(kRc);
    input_data_.b_cols = static_cast<int>(kRc);
    
    input_data_.a.assign(size, 0.0);
    input_data_.b.assign(size, 0.0);

    for (size_t i = 0; i < size; i++) {
        input_data_.a[i] = static_cast<double>(i % 100);
        input_data_.b[i] = static_cast<double>((i + 1) % 100);
    }

    expected_output_.assign(size, 0.0);

    for (size_t i = 0; i < kRc; ++i) {
        for (size_t k = 0; k < kRc; ++k) {
            double temp = input_data_.a[i * kRc + k];
            for (size_t j = 0; j < kRc; ++j) {
                expected_output_[i * kRc + j] += temp * input_data_.b[k * kRc + j];
            }
        }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }
    const double epsilon = 1e-7;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

private:
  InType input_data_;
  std::vector<double> expected_output_;
};

TEST_P(TabalaevAMatrixMulStrassenPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType,TabalaevAMatrixMulStrassenSEQ>(PPC_SETTINGS_tabalaev_a_matrix_mul_strassen);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TabalaevAMatrixMulStrassenPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TabalaevAMatrixMulStrassenPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tabalaev_a_matrix_mul_strassen
