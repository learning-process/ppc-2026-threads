#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <numeric>
#include <vector>

#include "smetanin_v_hoare_even_odd_batchelor/common/include/common.hpp"
#include "smetanin_v_hoare_even_odd_batchelor/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace smetanin_v_hoare_even_odd_batchelor {

class SmetaninVRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 1000000;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(static_cast<std::size_t>(kSize_));
    std::iota(input_data_.rbegin(), input_data_.rend(), 1);
  }

  void SetPerfAttributes(ppc::performance::PerfAttr &perf_attrs) override {
    const auto t0 = std::chrono::high_resolution_clock::now();
    perf_attrs.current_timer = [t0] {
      const auto now = std::chrono::high_resolution_clock::now();
      const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now - t0).count();
      return static_cast<double>(ns) * 1e-9;
    };
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SmetaninVRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SmetaninVHoarSortSEQ>(PPC_SETTINGS_smetanin_v_hoare_even_odd_batchelor);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SmetaninVRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SmetaninVRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace smetanin_v_hoare_even_odd_batchelor
