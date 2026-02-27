#include <gtest/gtest.h>

#include <cstddef>
#include <random>

#include "krymova_k_lsd_sort_merge_double/common/include/common.hpp"
#include "krymova_k_lsd_sort_merge_double/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace krymova_k_lsd_sort_merge_double {

class KrymovaKPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  const int size_ = 1000000;  // 1 миллион элементов для теста производительности

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10000.0, 10000.0);

    input_data_.resize(size_);
    for (int i = 0; i < size_; ++i) {
      input_data_[i] = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    for (size_t i = 1; i < output_data.size(); ++i) {
      if (output_data[i] < output_data[i - 1]) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KrymovaKPerfTests, MeasurePerformance) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, KrymovaKLsdSortMergeDoubleSEQ>(
    PPC_SETTINGS_krymova_k_lsd_sort_merge_double);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KrymovaKPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(PerformanceTests, KrymovaKPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace krymova_k_lsd_sort_merge_double