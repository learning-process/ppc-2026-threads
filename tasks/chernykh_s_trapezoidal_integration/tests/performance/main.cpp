#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <numbers>
#include <tuple>
#include <utility>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"
#include "chernykh_s_trapezoidal_integration/omp/include/ops_omp.hpp"
#include "chernykh_s_trapezoidal_integration/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chernykh_s_trapezoidal_integration {

class ChernykhSRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    limits_ = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
    steps_ = {400, 400, 400};
    auto f = [](const std::vector<double> &x) -> double { return std::sin(x[0]) * std::cos(x[1]) * std::exp(x[2]); };
    input_data_ = InType(limits_, steps_, f);
    reference_value_ = (1.0 - std::cos(1.0)) * std::sin(1.0) * (std::numbers::e - 1.0);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(reference_value_ - output_data) < 1e-1;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  std::vector<std::pair<double, double>> limits_;
  std::vector<int> steps_;
  InType input_data_{{}, {}, nullptr};
  OutType reference_value_ = 0.0;
};

TEST_P(ChernykhSRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, ChernykhSTrapezoidalIntegrationSEQ>(
                                              PPC_SETTINGS_chernykh_s_trapezoidal_integration),
                                          ppc::util::MakeAllPerfTasks<InType, ChernykhSTrapezoidalIntegrationOMP>(
                                              PPC_SETTINGS_chernykh_s_trapezoidal_integration));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChernykhSRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(IntegrationPerfTests, ChernykhSRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chernykh_s_trapezoidal_integration
