#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"
#include "zyuzin_n_multi_integrals_simpson/seq/include/ops_seq.hpp"

namespace zyuzin_n_multi_integrals_simpson {

class ZyuzinNRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    input_data_.lower_bounds = {0.0, 0.0, 0.0};
    input_data_.upper_bounds = {std::numbers::pi, std::numbers::pi / 2.0, 1.0};
    input_data_.n_steps = {180, 180, 120};
    input_data_.func = [](const std::vector<double> &p) { return std::sin(p[0]) * std::cos(p[1]) * std::exp(p[2]); };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double expected = 2.0 * (std::numbers::e - 1.0);
    return std::abs(output_data - expected) < 1e-3;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ZyuzinNRunPerfTestThreads, SimpsonTestRunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, ZyuzinNSimpsonSEQ>(PPC_SETTINGS_zyuzin_n_multi_integrals_simpson);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZyuzinNRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZyuzinNRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zyuzin_n_multi_integrals_simpson
