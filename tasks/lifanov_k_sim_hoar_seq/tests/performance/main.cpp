#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "lifanov_k_simple_hoar_seq/common/include/common.hpp"
#include "lifanov_k_simple_hoar_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace lifanov_k_simple_hoar_seq {

class LifanovKRunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const size_t kCount = 100000;
    input_data_.resize(kCount);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (size_t i = 0; i < kCount; ++i) {
      input_data_[i] = static_cast<int>(gen() % 10000);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::is_sorted(output_data.begin(), output_data.end());
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(LifanovKRunPerfTestSEQ, SortPerfMode) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, LifanovKSimpleHoarSEQ>("lifanov_k_simple_hoar_seq");

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LifanovKRunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(HoarSortPerf, LifanovKRunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace lifanov_k_simple_hoar_seq
