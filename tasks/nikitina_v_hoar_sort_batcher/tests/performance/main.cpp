#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"
#include "nikitina_v_hoar_sort_batcher/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nikitina_v_hoar_sort_batcher {

class NikitinaVHoarSortBatcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int count = 500000;
    input_data_.resize(count);

    std::mt19937 gen(42);
    std::uniform_int_distribution<> dist(-10000, 10000);
    for (int &x : input_data_) {
      x = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty() && std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(NikitinaVHoarSortBatcherPerfTests, RunPerfTests) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, HoareSortBatcherSEQ>(PPC_SETTINGS_nikitina_v_hoar_sort_batcher);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NikitinaVHoarSortBatcherPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(NikitinaVHoarSortBatcherPerfTests, NikitinaVHoarSortBatcherPerfTests, kGtestValues,
                         kPerfTestName);

}  // namespace
}  // namespace nikitina_v_hoar_sort_batcher
