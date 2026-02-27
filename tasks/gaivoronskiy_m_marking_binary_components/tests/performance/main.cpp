#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "gaivoronskiy_m_marking_binary_components/common/include/common.hpp"
#include "gaivoronskiy_m_marking_binary_components/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace gaivoronskiy_m_marking_binary_components {

class GaivoronskiyMMarkingPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int kSize = 500;
    input_data_.resize(static_cast<size_t>(kSize * kSize + 2));
    input_data_[0] = kSize;
    input_data_[1] = kSize;
    for (int i = 0; i < kSize; i++) {
      for (int j = 0; j < kSize; j++) {
        input_data_[static_cast<size_t>(i * kSize + j + 2)] = (i % 2 == 0) ? 0 : 1;
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    int rows = input_data_[0];
    int cols = input_data_[1];
    if (static_cast<int>(output_data.size()) != rows * cols + 2) {
      return false;
    }
    return std::any_of(output_data.begin() + 2, output_data.end(), [](int v) { return v > 0; });
  }

  InType GetTestInputData() final { return input_data_; }
};

TEST_P(GaivoronskiyMMarkingPerfTests, RunPerfModes) { ExecuteTest(GetParam()); }

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GaivoronskiyMMarkingBinaryComponentsSEQ>(
    PPC_SETTINGS_gaivoronskiy_m_marking_binary_components);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GaivoronskiyMMarkingPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GaivoronskiyMMarkingPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace gaivoronskiy_m_marking_binary_components
