#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>

#include "mityaeva_radix/common/include/common.hpp"
#include "mityaeva_radix/common/include/test_generator.hpp"
#include "mityaeva_radix/omp/include/ops_omp.hpp"
#include "mityaeva_radix/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace mityaeva_radix {

class MityaevaRadixPerf : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const std::size_t kCount_ = 4'000'000;
  InType input_data_;

  void SetUp() override {
    input_data_ = GenerateTest(kCount_, kCount_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::ranges::is_sorted(output_data, [](auto a, auto b) { return a <= b; });
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(MityaevaRadixPerf, RadixSortPerf) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MityaevaRadixSeq, MityaevaRadixOmp>(PPC_SETTINGS_mityaeva_radix);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MityaevaRadixPerf::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, MityaevaRadixPerf, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace mityaeva_radix
