#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <random>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"
#include "Nazarova_K_rad_sort_batcher_metod/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nazarova_k_rad_sort_batcher_metod_processes {

class NazarovaKRadSortBatcherMetodRunPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kCount = 200000;
  InType input_data_;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1e6, 1e6);

    input_data_.resize(static_cast<std::size_t>(kCount));
    for (int i = 0; i < kCount; ++i) {
      input_data_[static_cast<std::size_t>(i)] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    OutType expected = input_data_;
    std::ranges::sort(expected);
    if (output_data.size() != expected.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output_data.size(); ++i) {
      if (output_data[i] != expected[i]) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(NazarovaKRadSortBatcherMetodRunPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NazarovaKRadSortBatcherMetodSEQ>(PPC_SETTINGS_Nazarova_K_rad_sort_batcher_metod);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NazarovaKRadSortBatcherMetodRunPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NazarovaKRadSortBatcherMetodRunPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nazarova_k_rad_sort_batcher_metod_processes
