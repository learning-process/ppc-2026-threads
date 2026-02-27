#include <gtest/gtest.h>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"
#include "kichanova_k_lin_system_by_conjug_grad/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

class KichanovaKRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 50;
  InType input_data_{};

  void SetUp() override {
    LinSystemData data;
    data.n = kCount_;
    data.epsilon = 1e-8;
    data.A.assign(data.n * data.n, 0.0);
    for (int i = 0; i < data.n; ++i) {
        data.A[i * data.n + i] = 4.0;
        if (i > 0) data.A[i * data.n + (i - 1)] = -1.0;
        if (i < data.n - 1) data.A[i * data.n + (i + 1)] = -1.0;
    }
    data.b.resize(data.n);
    for (int i = 0; i < data.n; ++i) {
        data.b[i] = 1.0;
    }

    input_data_ = data;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != static_cast<size_t>(input_data_.n)) return false;
    
    double residual_norm = 0.0;
    for (int i = 0; i < input_data_.n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < input_data_.n; ++j) {
            sum += input_data_.A[i * input_data_.n + j] * output_data[j];
        }
        double diff = sum - input_data_.b[i];
        residual_norm += diff * diff;
    }
    residual_norm = std::sqrt(residual_norm);
    double tolerance = input_data_.epsilon * std::sqrt(static_cast<double>(input_data_.n));
    return residual_norm < tolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KichanovaKRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KichanovaKLinSystemByConjugGradSEQ>(PPC_SETTINGS_kichanova_k_lin_system_by_conjug_grad);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KichanovaKRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KichanovaKRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kichanova_k_lin_system_by_conjug_grad
