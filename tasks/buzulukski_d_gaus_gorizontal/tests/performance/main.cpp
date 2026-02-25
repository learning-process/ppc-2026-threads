#include <gtest/gtest.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"
#include "buzulukski_d_gaus_gorizontal/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace buzulukski_d_gaus_gorizontal {

class BuzulukskiDGausGorizontalPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  const int k_count = 1000;
  InType input_data{};

  void SetUp() override {
    input_data = k_count;
  }
  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == 100;
  }
  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(BuzulukskiDGausGorizontalPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BuzulukskiDGausGorizontalSEQ>(PPC_SETTINGS_buzulukski_d_gaus_gorizontal);
const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(RunModeTests, BuzulukskiDGausGorizontalPerfTests, kGtestValues,
                         BuzulukskiDGausGorizontalPerfTests::CustomPerfTestName);
}  // namespace

}  // namespace buzulukski_d_gaus_gorizontal
