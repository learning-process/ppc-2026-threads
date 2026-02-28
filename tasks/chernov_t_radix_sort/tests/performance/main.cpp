#include <gtest/gtest.h>

#include <algorithm>
#include <climits>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "chernov_t_radix_sort/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int size = 10'000'000;

    input_data_.resize(size);

    unsigned int cur_val = 69;

    for (auto &val : input_data_) {
      cur_val = (1664525 * cur_val) + 32433233;
      val = static_cast<int>(cur_val);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(ChernovTRadixSortPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, ChernovTRadixSortSEQ>(PPC_SETTINGS_chernov_t_radix_sort);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ChernovTRadixSortPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ChernovTRadixSortPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chernov_t_radix_sort
