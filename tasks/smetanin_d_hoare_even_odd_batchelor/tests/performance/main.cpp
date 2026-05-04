#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>

#include "performance/include/performance.hpp"
#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"
#include "smetanin_d_hoare_even_odd_batchelor/omp/include/ops_omp.hpp"
#include "smetanin_d_hoare_even_odd_batchelor/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace smetanin_d_hoare_even_odd_batchelor {

class SmetaninDRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 1000000;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(static_cast<std::size_t>(kSize_));
    for (std::size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = kSize_ - static_cast<int>(i);
    }
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
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SmetaninDRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SmetaninDHoarSortOMP, SmetaninDHoarSortSEQ>(
    PPC_SETTINGS_smetanin_d_hoare_even_odd_batchelor);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SmetaninDRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SmetaninDRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace smetanin_d_hoare_even_odd_batchelor
