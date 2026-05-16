#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"
#include "kruglova_a_conjugate_gradient_sle/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace kruglova_a_conjugate_gradient_sle {

class KruglovaAPerfTestAConjGradSle : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_size = 2000;
  InType input_data{};
  void SetUp() override {
    input_data.size = k_size;
    input_data.A.resize(static_cast<size_t>(k_size) * static_cast<size_t>(k_size));
    input_data.b.resize(k_size);

    for (int i = 0; i < k_size; ++i) {
      for (int j = 0; j < k_size; ++j) {
        if (i == j) {
          input_data.A[(i * k_size) + j] = static_cast<double>(k_size + 10);
        } else {
          input_data.A[(i * k_size) + j] = 1.0;
        }
      }
      input_data.b[i] = static_cast<double>((i % 10) + 1);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == static_cast<size_t>(k_size);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(KruglovaAPerfTestAConjGradSle, RunPerfPipeline) {
  ExecuteTest(GetParam());
}

TEST_P(KruglovaAPerfTestAConjGradSle, RunPerfTask) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KruglovaAConjGradSleSEQ

                                                       >(PPC_SETTINGS_kruglova_a_conjugate_gradient_sle);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KruglovaAPerfTestAConjGradSle::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SequentialPerformance, KruglovaAPerfTestAConjGradSle, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kruglova_a_conjugate_gradient_sle
