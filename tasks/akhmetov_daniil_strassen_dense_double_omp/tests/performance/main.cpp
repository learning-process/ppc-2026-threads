#include <gtest/gtest.h>

#include <cstddef>
#include <random>

#include "akhmetov_daniil_strassen_dense_double_omp/common/include/common.hpp"
#include "akhmetov_daniil_strassen_dense_double_omp/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace akhmetov_daniil_strassen_dense_double_omp {

class AkhmetovDaniilRunPerfTestsOMP : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const TestType n = 1024;

    input_data_.resize(1 + (2 * n * n));
    input_data_[0] = static_cast<double>(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (size_t i = 1; i < input_data_.size(); ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType& output_data) override {
    if (output_data.empty()) {
      return false;
    }

    const size_t n = format::GetN(input_data_);
    return output_data.size() == n * n;
  }

  InType GetTestInputData() override { return input_data_; }

 private:
  InType input_data_;
};

namespace {

TEST_P(AkhmetovDaniilRunPerfTestsOMP, StrassenTestPerf) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, AkhmetovDStrassenDenseDoubleOMP>(
    PPC_SETTINGS_akhmetov_daniil_strassen_dense_double_omp);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = AkhmetovDaniilRunPerfTestsOMP::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunStrassenPerfTestsOMP, AkhmetovDaniilRunPerfTestsOMP, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace akhmetov_daniil_strassen_dense_double_omp

