#include <gtest/gtest.h>

#include <algorithm>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"
#include "nikitina_v_hoar_sort_batcher/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nikitina_v_hoar_sort_batcher {

class NikitinaVHoarSortBatcherPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int count = 500000;
    input_data_.resize(count);

    int seed_val = 42;
    for (int &x : input_data_) {
      x = (seed_val % 20001) - 10000;
      seed_val = (seed_val * 73 + 17) % 20001;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty() && std::ranges::is_sorted(output_data);
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
