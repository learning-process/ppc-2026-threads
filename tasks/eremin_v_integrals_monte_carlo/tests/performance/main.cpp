#include <gtest/gtest.h>


#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"
#include "eremin_v_integrals_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace eremin_v_integrals_monte_carlo {

class EreminVRunPerfTestsThreadsIntegralsMonteCarlo : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
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

TEST_P(EreminVRunPerfTestsThreadsIntegralsMonteCarlo, IntegralsMonteCarloPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, EreminVIntegralsMonteCarloSEQ,
                                >(PPC_SETTINGS_eremin_v_strongin_algorithm);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = EreminVRunPerfTestsThreadsIntegralsMonteCarlo::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(IntegralsMonteCarloTestsPerf, EreminVRunPerfTestsThreadsIntegralsMonteCarlo, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace eremin_v_integrals_monte_carlo
