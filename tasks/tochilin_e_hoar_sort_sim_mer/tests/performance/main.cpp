#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"
#include "tochilin_e_hoar_sort_sim_mer/omp/include/ops_omp.hpp"
#include "tochilin_e_hoar_sort_sim_mer/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

class TochilinEHoarSortSimMerRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int k_count = 2000000;
  InType input_data;

  void SetUp() override {
    input_data.resize(static_cast<std::size_t>(k_count));
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dis(-10000, 10000);

    for (int i = 0; i < k_count; ++i) {
      input_data[static_cast<std::size_t>(i)] = dis(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data;
  }
};

TEST_P(TochilinEHoarSortSimMerRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, TochilinEHoarSortSimMerSEQ, TochilinEHoarSortSimMerOMP>(
    PPC_SETTINGS_tochilin_e_hoar_sort_sim_mer);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TochilinEHoarSortSimMerRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TochilinEHoarSortSimMerRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tochilin_e_hoar_sort_sim_mer
